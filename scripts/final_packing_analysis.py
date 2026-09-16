from __future__ import annotations
import sys, math, json
from pathlib import Path
from collections import defaultdict
import numpy as np, pandas as pd
import shapely
from shapely.geometry import box
from shapely.affinity import rotate
from shapely.ops import unary_union
from shapely import wkt
from scipy.stats import spearmanr
import statsmodels.api as sm

sys.path.insert(0,str(Path(__file__).resolve().parent))
import hardening_analysis as ha

OUT=Path(__file__).resolve().parents[1]/'results'

def edge_angles(poly):
    # A compact, reproducible candidate set from the longest roof edges.
    edges=[]
    rings=[poly.exterior]+list(poly.interiors)
    for ring in rings:
        cs=list(ring.coords)
        for a,b in zip(cs[:-1],cs[1:]):
            dx=b[0]-a[0];dy=b[1]-a[1];L=math.hypot(dx,dy)
            if L>=0.75:
                edges.append((L,math.degrees(math.atan2(dy,dx))%180.0))
    edges=sorted(edges,reverse=True)
    vals=[]
    for L,ang in edges:
        if all(min(abs(ang-v),180-abs(ang-v))>2 for v in vals):
            vals.append(ang)
        if len(vals)>=3: break
    vals.extend([ha.dominant_angle(poly)%180])
    out=[]
    for v in vals:
        out.extend([v,(v+90)%180])
    return sorted(set(round(v*2)/2 for v in out))

def count_angle(poly, angle, mod, steps=5):
    rp=rotate(poly,-angle,origin=poly.centroid,use_radians=False)
    minx,miny,maxx,maxy=rp.bounds; mw=mod['W'];ml=mod['L'];sx=mw+ha.GAP;sy=ml+ha.GAP
    best=0
    for ox in np.linspace(0,sx*.95,steps):
        xs=np.arange(minx+ox,maxx-mw+1e-9,sx)
        if xs.size==0: continue
        for oy in np.linspace(0,sy*.95,steps):
            ys=np.arange(miny+oy,maxy-ml+1e-9,sy)
            if ys.size==0: continue
            xx,yy=np.meshgrid(xs,ys); rects=shapely.box(xx.ravel(),yy.ravel(),xx.ravel()+mw,yy.ravel()+ml)
            n=int(np.count_nonzero(shapely.covers(rp,rects)))
            if n>best: best=n
    return best

def pack_edge(poly,mod,steps=4):
    if poly.is_empty:return 0
    parts=[poly] if poly.geom_type=='Polygon' else list(poly.geoms)
    total=0
    for p in parts:
        if p.area<mod['L']*mod['W']:continue
        best=0
        for ang in edge_angles(p):best=max(best,count_angle(p,ang,mod,steps))
        total+=best
    return total

def pack_exhaustive(poly,mod,angle_step=2,steps=7):
    if poly.is_empty:return 0
    parts=[poly] if poly.geom_type=='Polygon' else list(poly.geoms)
    total=0
    for p in parts:
        if p.area<mod['L']*mod['W']:continue
        best=0
        for ang in np.arange(0,180,angle_step):best=max(best,count_angle(p,float(ang),mod,steps))
        total+=best
    return total

def scenario(rows,setback,module_key):
    mod=ha.MODULES[module_key]; byb=defaultdict(list)
    for r in rows:byb[r['building_id']].append(r)
    rec=[]
    for bid,rs in byb.items():
        plans=[wkt.loads(r['plan_wkt']) for r in rs];env=unary_union(plans).buffer(0)
        shared=0
        for i in range(len(plans)):
            for j in range(i+1,len(plans)):
                shared+=plans[i].boundary.intersection(plans[j].boundary).length
        areas=np.array([r['area_3d'] for r in rs]); comps=[4*np.pi*p.area/p.length**2 if p.length else np.nan for p in plans]
        ua=0.;mods=0
        for r in rs:
            u=ha.usable_local_for_row(r,env,setback);ua+=u.area
            if not u.is_empty:mods+=pack_edge(u,mod,steps=4)
        area_kwp=ua*(mod['kWp']/(mod['L']*mod['W']));packed=mods*mod['kWp']
        small=areas[areas<2*mod['L']*mod['W']].sum()/areas.sum() if areas.sum() else np.nan
        rec.append({'building_id':bid,'building_name':rs[0]['building_name'],'n_facets':len(rs),'usable_area_m2':ua,'modules':mods,'packed_kwp':packed,'continuous_kwp':area_kwp,
                    'packing_shortfall_pct':100*(1-packed/area_kwp) if area_kwp else np.nan,'facets_per_100m2':100*len(rs)/ua if ua else np.nan,
                    'internal_boundary_density_m_per_m2':shared/env.area if env.area else np.nan,'shared_boundary_m':shared,'median_facet_area_m2':np.median(areas),
                    'facet_area_cv':areas.std()/areas.mean() if areas.mean() else np.nan,'small_facet_area_share':small,'mean_compactness':np.nanmean(comps),
                    'roof_envelope_area_m2':env.area,'setback_m':setback,'module':module_key})
    return pd.DataFrame(rec)

def main():
    raw=ha.parse_citygml();keep=[];seen=set()
    for r in raw:
        if (not r.get('repaired_valid',r['valid'])) or r['area_plan']<=.05:continue
        if (r['building_id'],r['dupkey']) in seen:continue
        seen.add((r['building_id'],r['dupkey']));keep.append(r)
    scenarios=[(0.0,ha.BASE_MODULE),(0.3,ha.BASE_MODULE),(0.5,ha.BASE_MODULE),(0.3,'LONGi_LR5_72HGD_575M')]
    summaries=[];base=None
    for sb,m in scenarios:
        print('final scenario',sb,m,flush=True);d=scenario(keep,sb,m);d.to_csv(OUT/f'FINAL_buildings_{m}_setback_{sb:.1f}.csv',index=False)
        if sb==.3 and m==ha.BASE_MODULE:base=d
        summaries.append({'setback_m':sb,'module':m,'n_buildings':len(d),'usable_area_m2':d.usable_area_m2.sum(),'continuous_kwp':d.continuous_kwp.sum(),'packed_kwp':d.packed_kwp.sum(),
                          'aggregate_shortfall_pct':100*(1-d.packed_kwp.sum()/d.continuous_kwp.sum()),'median_shortfall_pct':d.packing_shortfall_pct.median(),
                          'q25':d.packing_shortfall_pct.quantile(.25),'q75':d.packing_shortfall_pct.quantile(.75)})
    pd.DataFrame(summaries).to_csv(OUT/'FINAL_sensitivity_summary.csv',index=False)
    # correlations
    vars=['n_facets','facets_per_100m2','internal_boundary_density_m_per_m2','median_facet_area_m2','small_facet_area_share','facet_area_cv','mean_compactness','usable_area_m2']
    cr=[]
    for v in vars:
        rho,p=spearmanr(base[v],base.packing_shortfall_pct,nan_policy='omit');cr.append({'predictor':v,'rho':rho,'p':p})
    pd.DataFrame(cr).to_csv(OUT/'FINAL_correlations.csv',index=False)
    # regression HC3
    reg=base[['packing_shortfall_pct','internal_boundary_density_m_per_m2','usable_area_m2']].dropna().copy();reg['log_area']=np.log(reg.usable_area_m2)
    for c in ['internal_boundary_density_m_per_m2','log_area']:reg[c+'_z']=(reg[c]-reg[c].mean())/reg[c].std(ddof=0)
    model=sm.OLS(reg.packing_shortfall_pct,sm.add_constant(reg[['internal_boundary_density_m_per_m2_z','log_area_z']])).fit(cov_type='HC3')
    pd.DataFrame({'term':model.params.index,'coef':model.params.values,'se_HC3':model.bse.values,'p':model.pvalues.values}).to_csv(OUT/'FINAL_regression_HC3.csv',index=False)
    # leave-largest-out
    order=base.sort_values('usable_area_m2',ascending=False);loo=[]
    for k in [0,1,3]:
        d=order.iloc[k:];loo.append({'removed_largest_n':k,'n':len(d),'aggregate_shortfall_pct':100*(1-d.packed_kwp.sum()/d.continuous_kwp.sum())})
    pd.DataFrame(loo).to_csv(OUT/'FINAL_leave_largest_out.csv',index=False)
    # benchmark edge search against exhaustive on 10 moderate facets
    byb=defaultdict(list)
    for r in keep:byb[r['building_id']].append(r)
    cands=[]
    for bid,rs in byb.items():
        env=unary_union([wkt.loads(r['plan_wkt']) for r in rs]).buffer(0)
        for r in rs:
            u=ha.usable_local_for_row(r,env,.3)
            if 15<=u.area<=180:cands.append((u.area,bid,r['roof_id'],u))
    cands=sorted(cands,key=lambda x:x[0]);idx=np.linspace(0,len(cands)-1,10).round().astype(int);bench=[];mod=ha.MODULES[ha.BASE_MODULE]
    for i in idx:
        a,bid,rid,u=cands[i];edge=pack_edge(u,mod,steps=4);ref=pack_exhaustive(u,mod,angle_step=2,steps=7)
        bench.append({'building_id':bid,'roof_id':rid,'usable_area_m2':a,'edge_search':edge,'exhaustive_reference':ref,'gap_pct':100*(ref-edge)/ref if ref else 0});print('val',a,edge,ref,flush=True)
    bd=pd.DataFrame(bench);bd.to_csv(OUT/'FINAL_algorithm_validation.csv',index=False)
    summary={'retained_surfaces':len(keep),'buildings':len(base),'aggregate_shortfall_pct':100*(1-base.packed_kwp.sum()/base.continuous_kwp.sum()),'median_shortfall_pct':base.packing_shortfall_pct.median(),
             'q25':base.packing_shortfall_pct.quantile(.25),'q75':base.packing_shortfall_pct.quantile(.75),'usable_area_m2':base.usable_area_m2.sum(),'packed_kwp':base.packed_kwp.sum(),'continuous_kwp':base.continuous_kwp.sum(),
             'algorithm_validation_median_gap_pct':bd.gap_pct.median(),'algorithm_validation_max_gap_pct':bd.gap_pct.max()}
    (OUT/'FINAL_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':
    raise SystemExit('Internal implementation module. Use run_scenario.py, run_benchmark.py, or reproduce_empirical.py.')
