from __future__ import annotations
import sys, math, json, time
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np, pandas as pd
import shapely
from shapely import wkt, set_precision
from shapely.geometry import box
from shapely.affinity import rotate
from shapely.ops import unary_union
from scipy.stats import spearmanr
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
sys.path.insert(0,str(Path(__file__).resolve().parent))
import hardening_analysis as ha
import final_packing_analysis as fp

OUT=Path(__file__).resolve().parents[1]/'results'; OUT.mkdir(parents=True, exist_ok=True)
SEED=20260913
GRID_PRECISION=0.01
BASE_MODULE=ha.BASE_MODULE


def keep_rows():
    raw=ha.parse_citygml(); keep=[]; excluded=[]; seen=set()
    for r in raw:
        reason=None
        if not r.get('repaired_valid',r['valid']): reason='invalid_unrepairable'
        elif r['area_plan']<=0.05: reason='plan_area_le_0.05_m2'
        elif (r['building_id'],r['dupkey']) in seen: reason='exact_duplicate'
        if reason: excluded.append({**{k:v for k,v in r.items() if k not in ('pts3','basis')},'reason':reason})
        else:
            seen.add((r['building_id'],r['dupkey'])); keep.append(r)
    return raw, keep, excluded


def tilt_deg(row):
    basis=row.get('basis')
    if basis is None:return np.nan
    n=np.asarray(basis[3],float)
    if n[2]<0:n=-n
    return float(np.degrees(np.arccos(np.clip(n[2]/np.linalg.norm(n),-1,1))))


def edge_angles(poly):
    return fp.edge_angles(poly)


def count_angle_gap(poly, angle, mod, gap, steps):
    rp=rotate(poly,-angle,origin=poly.centroid,use_radians=False)
    minx,miny,maxx,maxy=rp.bounds; mw=mod['W']; ml=mod['L']; sx=mw+gap; sy=ml+gap
    best=0
    for ox in np.linspace(0,sx*.98,steps):
        xs=np.arange(minx+ox,maxx-mw+1e-9,sx)
        if xs.size==0: continue
        for oy in np.linspace(0,sy*.98,steps):
            ys=np.arange(miny+oy,maxy-ml+1e-9,sy)
            if ys.size==0: continue
            xx,yy=np.meshgrid(xs,ys)
            rects=shapely.box(xx.ravel(),yy.ravel(),xx.ravel()+mw,yy.ravel()+ml)
            n=int(np.count_nonzero(shapely.covers(rp,rects)))
            if n>best:best=n
    return best


def pack_production(poly,mod,gap=0.01):
    if poly.is_empty:return 0
    parts=[poly] if poly.geom_type=='Polygon' else list(poly.geoms)
    total=0
    for p in parts:
        if p.area<mod['L']*mod['W']:continue
        # Uniform search effort avoids shape-dependent search bias.
        steps=9
        angles=edge_angles(p)
        best=0
        for ang in angles:
            best=max(best,count_angle_gap(p,float(ang),mod,gap,steps))
        total+=best
    return total


def facet_worker(payload):
    bid,r,env_wkt,setback,module_key,gap=payload
    env=wkt.loads(env_wkt); mod=ha.MODULES[module_key]
    u=ha.usable_local_for_row(r,env,setback)
    return {'building_id':bid,'roof_id':r['roof_id'],'usable_area_m2':float(u.area),'modules':pack_production(u,mod,gap) if not u.is_empty else 0,
            'roof_area_3d_m2':r['area_3d'],'roof_area_plan_m2':r['area_plan'],'tilt_deg':tilt_deg(r)}


def scenario(keep,setback=.3,module_key=BASE_MODULE,gap=.02,exclude_flat_lt=None,workers=32):
    mod=ha.MODULES[module_key]; byb=defaultdict(list)
    for r in keep:byb[r['building_id']].append(r)
    tasks=[]; meta={}
    for bid,rs0 in byb.items():
        # flat-roof sensitivity removes low-tilt facets from BOTH packed and continuous reference.
        rs=[r for r in rs0 if exclude_flat_lt is None or tilt_deg(r)>=exclude_flat_lt]
        if not rs: continue
        plans0=[set_precision(wkt.loads(r['plan_wkt']),GRID_PRECISION) for r in rs]
        env=unary_union(plans0).buffer(0)
        total_perim=sum(p.length for p in plans0)
        shared=max(0.0,0.5*(total_perim-env.boundary.length))
        areas3d=np.array([r['area_3d'] for r in rs],float)
        comps=np.array([4*np.pi*p.area/p.length**2 if p.length else np.nan for p in plans0])
        gross3d=float(areas3d.sum()); grossplan=float(sum(p.area for p in plans0)); unionplan=float(env.area)
        small_thr=2*mod['L']*mod['W']; smallshare=float(areas3d[areas3d<small_thr].sum()/gross3d) if gross3d else np.nan
        meta[bid]={'building_name':rs[0]['building_name'],'n_facets':len(rs),'gross_roof_area_3d_m2':gross3d,'sum_plan_roof_area_m2':grossplan,'roof_envelope_area_m2':unionplan,
                   'plan_overlap_m2':grossplan-unionplan,'shared_boundary_m':shared,'internal_boundary_density_m_per_m2_grossplan':shared/unionplan if unionplan else np.nan,
                   'facet_density_per_100m2_gross3d':100*len(rs)/gross3d if gross3d else np.nan,'median_facet_area_m2':float(np.median(areas3d)),
                   'facet_area_cv':float(areas3d.std(ddof=0)/areas3d.mean()) if areas3d.mean() else np.nan,'small_facet_area_share':smallshare,
                   'mean_compactness':float(np.nanmean(comps)),'n_flat_lt5':sum(tilt_deg(r)<5 for r in rs0),'flat_area3d_lt5_m2':sum(r['area_3d'] for r in rs0 if tilt_deg(r)<5)}
        ew=env.wkt
        tasks.extend((bid,r,ew,setback,module_key,gap) for r in rs)
    facets=[]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        fut=[ex.submit(facet_worker,t) for t in tasks]
        for f in as_completed(fut): facets.append(f.result())
    fd=pd.DataFrame(facets); rows=[]
    for bid,g in fd.groupby('building_id'):
        ua=float(g.usable_area_m2.sum()); mods=int(g.modules.sum()); cont=ua*mod['kWp']/(mod['L']*mod['W']); packed=mods*mod['kWp']; m=meta[bid]
        rows.append({'building_id':bid,**m,'usable_area_m2':ua,'modules':mods,'packed_kwp':packed,'continuous_kwp':cont,
                     'packing_shortfall_pct':100*(1-packed/cont) if cont else np.nan,'setback_m':setback,'module':module_key,'gap_m':gap,'flat_exclusion_deg':exclude_flat_lt})
    return pd.DataFrame(rows).sort_values('building_id').reset_index(drop=True),fd


def deterministic_validation_sample(keep,n=12):
    # Candidate base-case facets. Selection is deterministic across a 4x3 grid of area and compactness strata.
    byb=defaultdict(list)
    for r in keep:byb[r['building_id']].append(r)
    rec=[]
    for bid,rs in byb.items():
        env=unary_union([wkt.loads(r['plan_wkt']) for r in rs]).buffer(0)
        for r in rs:
            u=ha.usable_local_for_row(r,env,.3)
            if u.is_empty or u.area<10 or u.area>220: continue
            comp=4*math.pi*u.area/u.length**2 if u.length else np.nan
            rec.append({'building_id':bid,'roof_id':r['roof_id'],'area':u.area,'compactness':comp,'geom':u})
    df=pd.DataFrame([{k:v for k,v in x.items() if k!='geom'} for x in rec])
    df['area_stratum']=pd.qcut(df['area'].rank(method='first'),4,labels=False)
    df['comp_stratum']=pd.qcut(df['compactness'].rank(method='first'),3,labels=False)
    chosen=[]
    for a in range(4):
        for c in range(3):
            sub=df[(df.area_stratum==a)&(df.comp_stratum==c)]
            if sub.empty:continue
            ta=sub.area.median(); tc=sub.compactness.median()
            score=((sub.area-ta)/(sub.area.std(ddof=0)+1e-9))**2+((sub.compactness-tc)/(sub.compactness.std(ddof=0)+1e-9))**2
            chosen.append(sub.loc[score.idxmin()])
    # should be 12; if any empty fill deterministic area quantiles not already used
    if len(chosen)<n:
        used={(x.building_id,x.roof_id) for x in chosen}
        rem=df[~df.apply(lambda r:(r.building_id,r.roof_id) in used,axis=1)].sort_values('area')
        idx=np.linspace(0,len(rem)-1,n-len(chosen)).round().astype(int)
        chosen.extend([rem.iloc[i] for i in idx])
    # reconnect geometries
    geommap={(x['building_id'],x['roof_id']):x['geom'] for x in rec}
    return [(r.building_id,r.roof_id,float(r.area),float(r.compactness),geommap[(r.building_id,r.roof_id)]) for r in chosen[:n]]


def benchmark_validation(keep):
    mod=ha.MODULES[BASE_MODULE]; rows=[]
    sample=deterministic_validation_sample(keep,12)
    for k,(bid,rid,area,comp,u) in enumerate(sample,1):
        prod=pack_production(u,mod,.01)
        # Stress-reference: 2-degree orientations plus exact edge angles, 21x21 offsets.
        ref=0
        angles=sorted(set([float(x) for x in np.arange(0,180,2)]+[float(x) for x in edge_angles(u)]))
        for ang in angles:
            ref=max(ref,count_angle_gap(u,ang,mod,.01,21))
        rows.append({'selection_index':k,'building_id':bid,'roof_id':rid,'usable_area_m2':area,'compactness':comp,'production_count':prod,'stress_reference_count':ref,
                     'gap_modules':ref-prod,'gap_pct':100*(ref-prod)/ref if ref else 0,'selection_rule':'4 area strata x 3 compactness strata; cell medoid'})
        print('benchmark',k,round(area,2),round(comp,3),prod,ref,flush=True)
    return pd.DataFrame(rows)


def correlation_table(base):
    vars=[
        ('n_facets','Raw facet count'),
        ('facet_density_per_100m2_gross3d','Facet density per 100 m2 gross 3D roof area'),
        ('median_facet_area_m2','Median facet area'),('small_facet_area_share','Small-facet area share'),('facet_area_cv','Facet-area coefficient of variation'),
        ('mean_compactness','Mean compactness'),('gross_roof_area_3d_m2','Gross 3D roof area')]
    rows=[]
    for v,label in vars:
        r,p=spearmanr(base[v],base.packing_shortfall_pct,nan_policy='omit');rows.append({'predictor':v,'label':label,'rho':r,'p_raw':p})
    df=pd.DataFrame(rows); df['p_fdr_bh']=multipletests(df.p_raw,method='fdr_bh')[1]
    return df


def regressions(base):
    v='facet_density_per_100m2_gross3d'; label='facet_density_gross3d'
    reg=base[['packing_shortfall_pct',v,'gross_roof_area_3d_m2']].dropna().copy()
    reg['log_area']=np.log(reg.gross_roof_area_3d_m2)
    for c in [v,'log_area']:
        reg[c+'_z']=(reg[c]-reg[c].mean())/reg[c].std(ddof=0)
    X=sm.add_constant(reg[[v+'_z','log_area_z']])
    m=sm.OLS(reg.packing_shortfall_pct,X).fit(cov_type='HC3'); ci=m.conf_int()
    return pd.DataFrame({'model':label,'term':m.params.index,'coef_pp':m.params.values,'se_HC3':m.bse.values,'p_value':m.pvalues.values,'ci_low':ci.iloc[:,0].values,'ci_high':ci.iloc[:,1].values})

def bootstrap_stability(base,nboot=10000):
    rng=np.random.default_rng(SEED); n=len(base)
    vals={'median':[],'aggregate':[],'rho_facet_density':[]}
    for _ in range(nboot):
        d=base.iloc[rng.integers(0,n,n)]
        vals['median'].append(d.packing_shortfall_pct.median())
        vals['aggregate'].append(100*(1-d.packed_kwp.sum()/d.continuous_kwp.sum()))
        vals['rho_facet_density'].append(spearmanr(d.facet_density_per_100m2_gross3d,d.packing_shortfall_pct,nan_policy='omit').statistic)
    est={'median':base.packing_shortfall_pct.median(),'aggregate':100*(1-base.packed_kwp.sum()/base.continuous_kwp.sum()),'rho_facet_density':spearmanr(base.facet_density_per_100m2_gross3d,base.packing_shortfall_pct).statistic}
    rows=[]
    for k,a in vals.items():
        rows.append({'metric':k,'estimate':est[k],'stability_low':np.nanquantile(a,.025),'stability_high':np.nanquantile(a,.975),'note':'within-window bootstrap stability interval; not a Stuttgart population confidence interval'})
    return pd.DataFrame(rows)

def mass_balance(keep,base):
    byb=defaultdict(list)
    for r in keep:byb[r['building_id']].append(r)
    sum_proj=sum(r['area_plan'] for r in keep); sum3d=sum(r['area_3d'] for r in keep)
    union_by_build=sum(unary_union([wkt.loads(r['plan_wkt']) for r in rs]).area for rs in byb.values())
    alkis=ha.parse_alkis(); alkis_union=unary_union([g for _,g in alkis]); alkis_sum=sum(g.area for _,g in alkis)
    return pd.DataFrame([
        {'quantity':'Export-log plan area','value_m2':71898.0,'note':'Metadata extent/area reported by supplied 2018 export log; not expected to equal roof-surface area'},
        {'quantity':'Sum of 58 AX_Gebaeude footprint areas','value_m2':alkis_sum,'note':'Cadastral features; overlaps/complexes can occur'},
        {'quantity':'Union of AX_Gebaeude footprints','value_m2':alkis_union.area,'note':'Plan-view cadastral footprint union'},
        {'quantity':'Sum of CityGML roof polygons in plan view','value_m2':sum_proj,'note':'Sum across 442 retained roof polygons; overlapping roof parts can inflate this relative to union'},
        {'quantity':'Sum of per-building CityGML roof-envelope unions','value_m2':union_by_build,'note':'Plan-view union within each of 31 CityGML buildings'},
        {'quantity':'Sum of retained CityGML 3D roof-surface areas','value_m2':sum3d,'note':'True surface area; can exceed plan area on pitched roofs'},
        {'quantity':'Base-case usable roof-plane area','value_m2':base.usable_area_m2.sum(),'note':'After 0.30 m planimetric external-envelope setback; used for continuous reference'}])


def main():
    raw,keep,ex=keep_rows(); print('raw/keep/ex',len(raw),len(keep),len(ex),flush=True)
    pd.DataFrame([{k:v for k,v in r.items() if k not in ('pts3','basis')} for r in ex]).to_csv(OUT/'geometry_exclusions.csv',index=False)
    scenarios=[('base',.3,BASE_MODULE,.01,None),('sb0',0,BASE_MODULE,.01,None),('sb05',.5,BASE_MODULE,.01,None),('large',.3,'LONGi_LR5_72HGD_575M',.01,None),
               ('gap0',.3,BASE_MODULE,0,None),('gap02',.3,BASE_MODULE,.02,None),('gap04',.3,BASE_MODULE,.04,None),('pitched_only5',.3,BASE_MODULE,.01,5.0)]
    summaries=[]; base=None
    for name,sb,mk,gap,flat in scenarios:
        t=time.time(); d,fd=scenario(keep,sb,mk,gap,flat,workers=32); d.to_csv(OUT/f'V5_{name}_buildings.csv',index=False); fd.to_csv(OUT/f'V5_{name}_facets.csv',index=False)
        if name=='base':base=d
        summaries.append({'scenario':name,'n_buildings':len(d),'setback_m':sb,'module':mk,'gap_m':gap,'exclude_flat_lt_deg':flat,'usable_area_m2':d.usable_area_m2.sum(),
                          'continuous_kwp':d.continuous_kwp.sum(),'packed_kwp':d.packed_kwp.sum(),'aggregate_shortfall_pct':100*(1-d.packed_kwp.sum()/d.continuous_kwp.sum()),
                          'median_building_shortfall_pct':d.packing_shortfall_pct.median(),'seconds':time.time()-t})
        print(name,summaries[-1]['aggregate_shortfall_pct'],summaries[-1]['seconds'],flush=True)
    pd.DataFrame(summaries).to_csv(OUT/'V5_sensitivity_summary.csv',index=False)
    cor=correlation_table(base); cor.to_csv(OUT/'V5_correlations_FDR.csv',index=False)
    reg=regressions(base); reg.to_csv(OUT/'V5_regressions_HC3.csv',index=False)
    boot=bootstrap_stability(base); boot.to_csv(OUT/'V5_bootstrap_stability.csv',index=False)
    # influence
    order=base.sort_values('gross_roof_area_3d_m2',ascending=False); loo=[]
    for k in [0,1,3]:
        d=order.iloc[k:]; loo.append({'removed_largest_n':k,'n_buildings':len(d),'aggregate_shortfall_pct':100*(1-d.packed_kwp.sum()/d.continuous_kwp.sum())})
    pd.DataFrame(loo).to_csv(OUT/'V5_leave_largest_out.csv',index=False)
    loor=[]
    for bid in base.building_id:
        d=base[base.building_id!=bid]
        loor.append({'removed_building_id':bid,'rho_facet_density':spearmanr(d.facet_density_per_100m2_gross3d,d.packing_shortfall_pct,nan_policy='omit').statistic})
    pd.DataFrame(loor).to_csv(OUT/'V5_leave_one_building_out_rho.csv',index=False)
    # existing alignment sensitivity from V3/consistent geometry; recompute merge with base
    al=pd.read_csv(str(Path(__file__).resolve().parents[1]/'results'/'alkis_citygml_alignment_bestsubset.csv'))
    # detect column
    ioucol=[c for c in al.columns if 'iou' in c.lower()][0]
    merged=base.merge(al[['building_id',ioucol]],on='building_id',how='left')
    arows=[]
    for th in [.8,.9]:
        d=merged[merged[ioucol]>=th]; arows.append({'iou_threshold':th,'n_buildings':len(d),'aggregate_shortfall_pct':100*(1-d.packed_kwp.sum()/d.continuous_kwp.sum())})
    pd.DataFrame(arows).to_csv(OUT/'V5_alignment_sensitivity.csv',index=False)
    # validation benchmark
    bench=benchmark_validation(keep); bench.to_csv(OUT/'V5_packing_stress_benchmark.csv',index=False)
    # mass balance and counts
    mb=mass_balance(keep,base); mb.to_csv(OUT/'V5_mass_balance.csv',index=False)
    counts=pd.DataFrame([
        {'stage':'Stuttgart export summary metadata','count':94,'unit':'Gebäude records','note':'Metadata count from export summary; not a sequential filtering stage'},
        {'stage':'Parsed NAS feature type AX_Gebaeude','count':58,'unit':'features','note':'Cadastral building-footprint objects in supplied NAS'},
        {'stage':'CityGML Building objects with roof geometry','count':31,'unit':'buildings','note':'Separate 3D model object count in supplied CityGML; not expected to equal NAS count'},
        {'stage':'Raw CityGML roof polygons','count':450,'unit':'polygons','note':'Before geometry audit'},
        {'stage':'Retained CityGML roof polygons','count':442,'unit':'polygons','note':'After repair/exclusion audit'}])
    counts.to_csv(OUT/'V5_data_lineage_counts.csv',index=False)
    summary={'base':summaries[0],'correlations':cor.to_dict('records'),'bootstrap_stability':boot.to_dict('records'),'validation_max_gap_modules':int(bench.gap_modules.max()),
             'validation_max_gap_pct':float(bench.gap_pct.max()),'flat_surface_count':int((pd.read_csv(OUT/'V5_base_facets.csv').tilt_deg<5).sum()),
             'flat_surface_area3d_m2':float(pd.read_csv(OUT/'V5_base_facets.csv').loc[lambda x:x.tilt_deg<5,'roof_area_3d_m2'].sum())}
    (OUT/'V5_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    raise SystemExit('Internal implementation module. Use run_scenario.py, run_benchmark.py, or reproduce_empirical.py.')
