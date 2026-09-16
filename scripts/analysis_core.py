from __future__ import annotations
import sys, math, json, hashlib, time
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np, pandas as pd
import shapely
from shapely import wkt, set_precision
from shapely.ops import unary_union
from shapely.affinity import rotate
from scipy.stats import spearmanr
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

SRC=Path(__file__).resolve().parent
sys.path.insert(0,str(SRC))
import hardening_analysis as ha
import final_revision_v5 as old

OUT=Path(__file__).resolve().parents[1]/'results'; OUT.mkdir(parents=True,exist_ok=True)
SEED=20260913
BASE=ha.BASE_MODULE
GRID=.01

# ---------- geometry cleaning with row-stable identifiers ----------
def plane_from_pts(pts):
    p=np.asarray(pts,float); c=p.mean(axis=0); _,_,vh=np.linalg.svd(p-c,full_matrices=False); n=vh[-1]
    if n[2]<0:n=-n
    n=n/np.linalg.norm(n); return n,-float(np.dot(n,c))

def z_at(x,y,n,d):
    return -(n[0]*x+n[1]*y+d)/n[2] if abs(n[2])>1e-10 else np.nan

def initial_rows():
    raw=ha.parse_citygml(); keep=[]; excluded=[]; seen=set()
    for idx,r0 in enumerate(raw):
        r=dict(r0); r['surface_uid']=f"S{idx+1:04d}"
        reason=None
        if not r.get('repaired_valid',r['valid']): reason='invalid_unrepairable'
        elif r['area_plan']<=0.05: reason='plan_area_le_0.05_m2'
        elif (r['building_id'],r['dupkey']) in seen: reason='exact_duplicate_geometry'
        if reason: excluded.append({'surface_uid':r['surface_uid'],'building_id':r['building_id'],'original_roof_id':r['roof_id'],'reason':reason})
        else:
            seen.add((r['building_id'],r['dupkey'])); keep.append(r)
    return raw,keep,excluded

def remove_nested(rows, angle_tol=1.0,z_tol=.10,contain=.99):
    by=defaultdict(list)
    for r in rows:by[r['building_id']].append(r)
    remove=set(); audit=[]
    for bid,rs in by.items():
        plans=[wkt.loads(r['plan_wkt']) for r in rs]; planes=[plane_from_pts(r['pts3']) for r in rs]
        for i in range(len(rs)):
            for j in range(i+1,len(rs)):
                pi,pj=plans[i],plans[j]; inter=pi.intersection(pj); a=inter.area
                if a<=.01:continue
                ni,di=planes[i]; nj,dj=planes[j]
                ang=float(np.degrees(np.arccos(np.clip(abs(float(np.dot(ni,nj))),-1,1))))
                rp=inter.representative_point(); zi=z_at(rp.x,rp.y,ni,di);zj=z_at(rp.x,rp.y,nj,dj);dz=abs(zi-zj) if np.isfinite(zi) and np.isfinite(zj) else np.nan
                frac=a/min(pi.area,pj.area)
                flag=ang<=angle_tol and dz<=z_tol and frac>=contain
                removed=''
                if flag:
                    k=i if pi.area<=pj.area else j; remove.add(rs[k]['surface_uid']); removed=rs[k]['surface_uid']
                audit.append({'building_id':bid,'uid_a':rs[i]['surface_uid'],'uid_b':rs[j]['surface_uid'],'overlap_m2':a,'overlap_fraction_smaller':frac,'angle_deg':ang,'vertical_separation_m':dz,'nested_coplanar':flag,'removed_uid':removed})
    kept=[r for r in rows if r['surface_uid'] not in remove]
    return kept,remove,pd.DataFrame(audit)

def merge_coplanar_adjacent(rows, angle_tol=1.0,z_tol=.10,shared_min=.05):
    """Merge contiguous polygons that lie on the same physical plane.
    The rule is purely geometric: >5 cm shared plan-view boundary, plane-normal difference <=1 deg,
    and vertical difference at the shared boundary <=0.10 m.
    """
    by=defaultdict(list)
    for r in rows:by[r['building_id']].append(r)
    merged=[]; pair_audit=[]; membership=[]
    for bid,rs in by.items():
        n=len(rs); parent=list(range(n))
        def find(a):
            while parent[a]!=a:
                parent[a]=parent[parent[a]]; a=parent[a]
            return a
        def union(a,b):
            ra,rb=find(a),find(b)
            if ra!=rb: parent[rb]=ra
        plans=[set_precision(wkt.loads(r['plan_wkt']),GRID) for r in rs]; planes=[plane_from_pts(r['pts3']) for r in rs]
        for i in range(n):
            for j in range(i+1,n):
                shared=plans[i].boundary.intersection(plans[j].boundary)
                L=shared.length
                if L<=shared_min:continue
                ni,di=planes[i];nj,dj=planes[j]
                ang=float(np.degrees(np.arccos(np.clip(abs(float(np.dot(ni,nj))),-1,1))))
                rp=shared.representative_point(); zi=z_at(rp.x,rp.y,ni,di);zj=z_at(rp.x,rp.y,nj,dj);dz=abs(zi-zj) if np.isfinite(zi) and np.isfinite(zj) else np.nan
                flag=ang<=angle_tol and dz<=z_tol
                pair_audit.append({'building_id':bid,'uid_a':rs[i]['surface_uid'],'uid_b':rs[j]['surface_uid'],'shared_boundary_m':L,'angle_deg':ang,'vertical_separation_m':dz,'merged_same_plane':flag,'same_building_part':rs[i]['part_id']==rs[j]['part_id']})
                if flag:union(i,j)
        groups=defaultdict(list)
        for i in range(n):groups[find(i)].append(i)
        for k,inds in enumerate(groups.values(),1):
            members=[rs[i]['surface_uid'] for i in inds]
            gplan=unary_union([plans[i] for i in inds]).buffer(0)
            # same-plane component: first basis valid for the full merged polygon
            basis=rs[inds[0]]['basis']; glocal=ha.plan_geom_to_local(gplan,basis)
            normal=np.asarray(basis[3],float); normal=normal/np.linalg.norm(normal)
            tilt=float(np.degrees(np.arccos(np.clip(abs(normal[2]),-1,1))))
            area3d=float(glocal.area)
            uid=f"{bid}_PF{k:03d}"
            rr={'building_id':bid,'building_name':rs[inds[0]]['building_name'],'physical_facet_uid':uid,
                'source_surface_uids':'|'.join(members),'source_surface_count':len(members),'part_ids':'|'.join(sorted(set(rs[i]['part_id'] or '' for i in inds))),
                'plan_wkt':gplan.wkt,'basis':basis,'area_plan':float(gplan.area),'area_3d':area3d,'tilt_deg':tilt}
            merged.append(rr)
            for m in members:membership.append({'building_id':bid,'surface_uid':m,'physical_facet_uid':uid,'component_size':len(members)})
    return merged,pd.DataFrame(pair_audit),pd.DataFrame(membership)

# ---------- packing ----------
def edge_angles(poly):return old.edge_angles(poly)
def count_angle_gap(poly,ang,mod,gap,steps):return old.count_angle_gap(poly,ang,mod,gap,steps)
def pack_production(poly,mod,gap=.01):return old.pack_production(poly,mod,gap)

def local_for_facet(r,env,setback):
    e=env.buffer(-setback,join_style=2) if setback>0 else env
    clipped=wkt.loads(r['plan_wkt']).intersection(e)
    return ha.plan_geom_to_local(clipped,r['basis'])

def facet_task(t):
    bid,r,envw,sb,mkey,gap=t; env=wkt.loads(envw);mod=ha.MODULES[mkey]
    u=local_for_facet(r,env,sb)
    return {'building_id':bid,'physical_facet_uid':r['physical_facet_uid'],'usable_area_m2':float(u.area),'modules':pack_production(u,mod,gap) if not u.is_empty else 0,
            'gross_area_3d_m2':r['area_3d'],'gross_area_plan_m2':r['area_plan'],'tilt_deg':r['tilt_deg'],'source_surface_count':r['source_surface_count']}

def scenario(facets,setback=.3,mkey=BASE,gap=.01,tilt_subset=None,workers=24):
    mod=ha.MODULES[mkey];by=defaultdict(list)
    for r in facets:by[r['building_id']].append(r)
    tasks=[];meta={}
    for bid,rs0 in by.items():
        if tilt_subset=='pitched':rs=[r for r in rs0 if r['tilt_deg']>=5]
        elif tilt_subset=='flat':rs=[r for r in rs0 if r['tilt_deg']<5]
        else:rs=rs0
        if not rs:continue
        plans=[set_precision(wkt.loads(r['plan_wkt']),GRID) for r in rs];env=unary_union(plans).buffer(0); env=wkt.loads(env.wkt)  # match production serialization and remove precision-model side effects
        areas=np.array([r['area_3d'] for r in rs],float)
        # Compactness is measured on the roof plane (not the horizontal projection).
        local_full=[ha.plan_geom_to_local(p,r['basis']) for p,r in zip(plans,rs)]
        comps=np.array([4*np.pi*g.area/(g.length**2) if g.length else np.nan for g in local_full])
        small_thr=2*mod['L']*mod['W']
        gross=float(areas.sum())
        meta[bid]={'building_name':rs[0]['building_name'],'n_facets':len(rs),'gross_roof_area_3d_m2':gross,
                   'facet_density_per_100m2_gross3d':100*len(rs)/gross if gross else np.nan,
                   'median_facet_area_m2':float(np.median(areas)),'facet_area_cv':float(areas.std(ddof=0)/areas.mean()) if areas.mean() else np.nan,
                   'small_facet_area_share':float(areas[areas<small_thr].sum()/gross) if gross else np.nan,
                   'mean_compactness':float(np.nanmean(comps)),'small_facet_threshold_m2':small_thr}
        ew=env.wkt;tasks += [(bid,r,ew,setback,mkey,gap) for r in rs]
    rec=[]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        fs=[ex.submit(facet_task,t) for t in tasks]
        for f in as_completed(fs):rec.append(f.result())
    fd=pd.DataFrame(rec);rows=[]
    for bid,g in fd.groupby('building_id'):
        ua=float(g.usable_area_m2.sum());mods=int(g.modules.sum());cont=ua*mod['kWp']/(mod['L']*mod['W']);packed=mods*mod['kWp']
        rows.append({'building_id':bid,**meta[bid],'usable_area_m2':ua,'modules':mods,'packed_kwp':packed,'continuous_kwp':cont,'packing_shortfall_pct':100*(1-packed/cont) if cont else np.nan})
    return pd.DataFrame(rows).sort_values('building_id').reset_index(drop=True),fd.sort_values(['building_id','physical_facet_uid']).reset_index(drop=True)

# ---------- statistics ----------
def corr_table(base):
    vs=[('n_facets','Raw facet count'),('facet_density_per_100m2_gross3d','Facet density per 100 m² gross 3D roof area'),('median_facet_area_m2','Median facet area'),('small_facet_area_share','Small-facet area share'),('facet_area_cv','Facet-area coefficient of variation'),('mean_compactness','Mean compactness'),('gross_roof_area_3d_m2','Gross 3D roof area')]
    rows=[]
    for c,label in vs:
        rho,p=spearmanr(base[c],base.packing_shortfall_pct,nan_policy='omit');rows.append({'predictor':c,'label':label,'rho':rho,'p_raw':p})
    d=pd.DataFrame(rows);d['p_fdr_bh']=multipletests(d.p_raw,method='fdr_bh')[1];return d

def regression(base):
    d=base[['packing_shortfall_pct','facet_density_per_100m2_gross3d','gross_roof_area_3d_m2']].dropna().copy();d['log_area']=np.log(d.gross_roof_area_3d_m2)
    for c in ['facet_density_per_100m2_gross3d','log_area']:d[c+'_z']=(d[c]-d[c].mean())/d[c].std(ddof=0)
    X=sm.add_constant(d[['facet_density_per_100m2_gross3d_z','log_area_z']]);m=sm.OLS(d.packing_shortfall_pct,X).fit(cov_type='HC3');ci=m.conf_int()
    return pd.DataFrame({'term':m.params.index,'coef_pp':m.params.values,'se_HC3':m.bse.values,'p_value':m.pvalues.values,'ci_low':ci.iloc[:,0].values,'ci_high':ci.iloc[:,1].values})

def bootstrap(base,n=10000):
    rng=np.random.default_rng(SEED);N=len(base);a=[];med=[];rho=[]
    for _ in range(n):
        d=base.iloc[rng.integers(0,N,N)];a.append(100*(1-d.packed_kwp.sum()/d.continuous_kwp.sum()));med.append(d.packing_shortfall_pct.median());rho.append(spearmanr(d.facet_density_per_100m2_gross3d,d.packing_shortfall_pct).statistic)
    rows=[]
    for name,x,est in [('aggregate',a,100*(1-base.packed_kwp.sum()/base.continuous_kwp.sum())),('median',med,base.packing_shortfall_pct.median()),('rho_facet_density',rho,spearmanr(base.facet_density_per_100m2_gross3d,base.packing_shortfall_pct).statistic)]:
        rows.append({'metric':name,'estimate':est,'percentile_2.5':np.nanquantile(x,.025),'percentile_97.5':np.nanquantile(x,.975),'n_resamples':n,'seed':SEED})
    return pd.DataFrame(rows)

# ---------- deterministic full-range benchmark with unique IDs ----------
def benchmark(facets):
    mod=ha.MODULES[BASE];rec=[];geom={}
    by=defaultdict(list)
    for r in facets:by[r['building_id']].append(r)
    for bid,rs in by.items():
        env=unary_union([wkt.loads(r['plan_wkt']) for r in rs]).buffer(0)
        for r in rs:
            u=local_for_facet(r,env,.3)
            if u.is_empty or u.area<mod['L']*mod['W']:continue
            comp=4*np.pi*u.area/u.length**2 if u.length else np.nan;uid=r['physical_facet_uid']
            rec.append({'building_id':bid,'surface_uid':uid,'area':u.area,'compactness':comp});geom[uid]=u
    d=pd.DataFrame(rec);d['area_stratum']=pd.qcut(d.area.rank(method='first'),5,labels=False);d['comp_stratum']=pd.qcut(d.compactness.rank(method='first'),3,labels=False)
    sel=[]
    for a in range(5):
        for c in range(3):
            s=d[(d.area_stratum==a)&(d.comp_stratum==c)].copy()
            if s.empty:continue
            za=(s.area-s.area.median())/(s.area.std(ddof=0)+1e-12);zc=(s.compactness-s.compactness.median())/(s.compactness.std(ddof=0)+1e-12)
            sel.append(s.loc[(za*za+zc*zc).idxmin()])
    out=[]
    for k,r in enumerate(sel,1):
        u=geom[r.surface_uid];prod=pack_production(u,mod,.01);ref=0
        # One benchmark definition only: complete 5° sweep plus exact edge-derived angles, 13x13 offsets.
        angs=sorted(set([float(x) for x in np.arange(0,180,5)]+[float(x) for x in edge_angles(u)]))
        for ang in angs:ref=max(ref,count_angle_gap(u,ang,mod,.01,13))
        out.append({'selection_index':k,'building_id':r.building_id,'surface_uid':r.surface_uid,'usable_area_m2':float(r.area),'compactness':float(r.compactness),'area_stratum':int(r.area_stratum),'compactness_stratum':int(r.comp_stratum),'production_count':prod,'stress_count':ref,'gap_modules':ref-prod,'gap_pct':100*(ref-prod)/ref if ref else 0,'reference_search':'5-degree sweep + exact production edge angles; 13x13 offset grid'})
    return pd.DataFrame(out)

def main():
    raw,k0,ex0=initial_rows();k1,removed,nest=remove_nested(k0);facets,pairs,membership=merge_coplanar_adjacent(k1)
    print('raw',len(raw),'after basic',len(k0),'after nested',len(k1),'physical facets',len(facets),flush=True)
    pd.DataFrame(ex0).to_csv(OUT/'geometry_basic_exclusions.csv',index=False);nest.to_csv(OUT/'nested_overlap_audit.csv',index=False);pairs.to_csv(OUT/'coplanar_adjacency_audit.csv',index=False);membership.to_csv(OUT/'surface_to_physical_facet.csv',index=False)
    pd.DataFrame([{k:v for k,v in r.items() if k!='basis'} for r in facets]).to_csv(OUT/'physical_facets.csv',index=False)
    scenarios=[('base',.3,BASE,.01,None),('setback0',0,BASE,.01,None),('setback05',.5,BASE,.01,None),('large575',.3,'LONGi_LR5_72HGD_575M',.01,None),('gap0',.3,BASE,0,None),('gap02',.3,BASE,.02,None),('gap04',.3,BASE,.04,None),('pitched',.3,BASE,.01,'pitched'),('flat',.3,BASE,.01,'flat')]
    sums=[];base=None
    for name,sb,m,g,t in scenarios:
        t0=time.time();bd,fd=scenario(facets,sb,m,g,t,workers=5);bd.to_csv(OUT/f'{name}_buildings.csv',index=False);fd.to_csv(OUT/f'{name}_facets.csv',index=False)
        if name=='base':base=bd
        sums.append({'scenario':name,'n_buildings':len(bd),'setback_m':sb,'module':m,'gap_m':g,'tilt_subset':t or 'all','usable_area_m2':bd.usable_area_m2.sum(),'continuous_kwp':bd.continuous_kwp.sum(),'packed_kwp':bd.packed_kwp.sum(),'aggregate_shortfall_pct':100*(1-bd.packed_kwp.sum()/bd.continuous_kwp.sum()),'median_building_shortfall_pct':bd.packing_shortfall_pct.median(),'seconds':time.time()-t0});print(name,sums[-1]['aggregate_shortfall_pct'],flush=True)
    pd.DataFrame(sums).to_csv(OUT/'sensitivity_summary.csv',index=False)
    cor=corr_table(base);cor.to_csv(OUT/'correlations_FDR.csv',index=False);reg=regression(base);reg.to_csv(OUT/'regression_HC3.csv',index=False);boot=bootstrap(base);boot.to_csv(OUT/'bootstrap_percentile.csv',index=False)
    # influence
    order=base.sort_values('gross_roof_area_3d_m2',ascending=False);pd.DataFrame([{'removed_largest_n':k,'aggregate_shortfall_pct':100*(1-order.iloc[k:].packed_kwp.sum()/order.iloc[k:].continuous_kwp.sum())} for k in [0,1,3]]).to_csv(OUT/'leave_largest_out.csv',index=False)
    loo=[]
    for bid in base.building_id:
        d=base[base.building_id!=bid];loo.append({'removed_building_id':bid,'rho':spearmanr(d.facet_density_per_100m2_gross3d,d.packing_shortfall_pct).statistic})
    pd.DataFrame(loo).to_csv(OUT/'leave_one_out.csv',index=False)
    # alignment sensitivity reuse exact building IoUs from prior deterministic check
    al=pd.read_csv(str(Path(__file__).resolve().parents[1]/'results'/'alkis_citygml_alignment_bestsubset.csv'));iou=[c for c in al if 'iou' in c.lower()][0];m=base.merge(al[['building_id',iou]],on='building_id',how='left');arr=[]
    for th in [.8,.9]:
        d=m[m[iou]>=th];arr.append({'iou_threshold':th,'n_buildings':len(d),'aggregate_shortfall_pct':100*(1-d.packed_kwp.sum()/d.continuous_kwp.sum())})
    pd.DataFrame(arr).to_csv(OUT/'alignment_sensitivity.csv',index=False)
    b=benchmark(facets);b.to_csv(OUT/'packing_benchmark.csv',index=False)
    # compact summary
    sumj={'raw_polygons':len(raw),'after_basic_cleaning':len(k0),'after_nested_removal':len(k1),'physical_facets_after_coplanar_merge':len(facets),'merged_adjacency_pairs':int(pairs.merged_same_plane.sum()),'buildings_affected_by_merge':int(membership.groupby('building_id').apply(lambda x:(x.component_size>1).any(),include_groups=False).sum()),'base':sums[0],'bootstrap':boot.to_dict('records'),'correlations':cor.to_dict('records'),'benchmark_agree':int((b.gap_modules==0).sum()),'benchmark_n':len(b),'benchmark_max_gap':int(b.gap_modules.max()),'benchmark_pooled_pct':100*(b.stress_count.sum()-b.production_count.sum())/b.stress_count.sum()}
    Path(OUT/'summary.json').write_text(json.dumps(sumj,indent=2));print(json.dumps(sumj,indent=2),flush=True)

if __name__=='__main__':
    raise SystemExit('Internal implementation module. Use run_scenario.py, run_benchmark.py, or reproduce_empirical.py.')

def angle_worker(payload):
    key,poly_wkt,ang,mkey,gap=payload
    p=wkt.loads(poly_wkt);mod=ha.MODULES[mkey]
    return key,count_angle_gap(p,ang,mod,gap,9)

def scenario_fast(facets,setback=.3,mkey=BASE,gap=.01,tilt_subset=None,workers=5):
    mod=ha.MODULES[mkey];by=defaultdict(list)
    for r in facets:by[r['building_id']].append(r)
    meta={};facet_info=[];tasks=[]
    for bid,rs0 in by.items():
        if tilt_subset=='pitched':rs=[r for r in rs0 if r['tilt_deg']>=5]
        elif tilt_subset=='flat':rs=[r for r in rs0 if r['tilt_deg']<5]
        else:rs=rs0
        if not rs:continue
        plans=[set_precision(wkt.loads(r['plan_wkt']),GRID) for r in rs];env=unary_union(plans).buffer(0); env=wkt.loads(env.wkt)  # match production serialization and remove precision-model side effects
        areas=np.array([r['area_3d'] for r in rs],float);local_full=[ha.plan_geom_to_local(p,r['basis']) for p,r in zip(plans,rs)]
        comps=np.array([4*np.pi*g.area/(g.length**2) if g.length else np.nan for g in local_full]);gross=float(areas.sum());small_thr=2*mod['L']*mod['W']
        meta[bid]={'building_name':rs[0]['building_name'],'n_facets':len(rs),'gross_roof_area_3d_m2':gross,'facet_density_per_100m2_gross3d':100*len(rs)/gross if gross else np.nan,'median_facet_area_m2':float(np.median(areas)),'facet_area_cv':float(areas.std(ddof=0)/areas.mean()) if areas.mean() else np.nan,'small_facet_area_share':float(areas[areas<small_thr].sum()/gross) if gross else np.nan,'mean_compactness':float(np.nanmean(comps)),'small_facet_threshold_m2':small_thr}
        e=env.buffer(-setback,join_style=2) if setback>0 else env
        for r in rs:
            clipped=wkt.loads(r['plan_wkt']).intersection(e);u=ha.plan_geom_to_local(clipped,r['basis']);ua=float(u.area)
            parts=[] if u.is_empty else ([u] if u.geom_type=='Polygon' else [g for g in u.geoms if g.geom_type=='Polygon'])
            part_keys=[]
            for pi,p in enumerate(parts):
                key=(r['physical_facet_uid'],pi);part_keys.append(key)
                if p.area>=mod['L']*mod['W']:
                    for ang in edge_angles(p):tasks.append((key,p.wkt,float(ang),mkey,gap))
            facet_info.append({'building_id':bid,'physical_facet_uid':r['physical_facet_uid'],'usable_area_m2':ua,'gross_area_3d_m2':r['area_3d'],'gross_area_plan_m2':r['area_plan'],'tilt_deg':r['tilt_deg'],'source_surface_count':r['source_surface_count'],'part_keys':part_keys})
    best=defaultdict(int)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        fs=[ex.submit(angle_worker,t) for t in tasks]
        for f in as_completed(fs):
            key,n=f.result();best[key]=max(best[key],n)
    rec=[]
    for f in facet_info:
        mods=sum(best[k] for k in f.pop('part_keys'));rec.append({**f,'modules':mods})
    fd=pd.DataFrame(rec);rows=[]
    for bid,g in fd.groupby('building_id'):
        ua=float(g.usable_area_m2.sum());mods=int(g.modules.sum());cont=ua*mod['kWp']/(mod['L']*mod['W']);packed=mods*mod['kWp']
        rows.append({'building_id':bid,**meta[bid],'usable_area_m2':ua,'modules':mods,'packed_kwp':packed,'continuous_kwp':cont,'packing_shortfall_pct':100*(1-packed/cont) if cont else np.nan})
    return pd.DataFrame(rows).sort_values('building_id').reset_index(drop=True),fd.sort_values(['building_id','physical_facet_uid']).reset_index(drop=True)
