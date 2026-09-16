from __future__ import annotations

import json, math, zipfile, hashlib
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from lxml import etree
from shapely.geometry import Polygon, MultiPolygon, box, Point
from shapely.ops import unary_union
from shapely.affinity import rotate
from shapely.prepared import prep
from shapely import wkt
import shapely
from scipy.stats import spearmanr
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'; OUT.mkdir(parents=True, exist_ok=True)
CITYGML_ZIP = ROOT / 'data_raw' / 'CityGML_Alles_2018_11_27_07_15_56.zip'
ALKIS_XML = ROOT / 'data_raw' / '2018-11-27_07-16-14_Liegenschaftskarte_ohne_Eigentuemer.xml'

MODULES = {
    'LONGi_LR5_54HTB_440M': {'L':1.722,'W':1.134,'kWp':0.440},
    'LONGi_LR5_72HGD_575M': {'L':2.278,'W':1.134,'kWp':0.575},
}
GAP=0.02
BASE_MODULE='LONGi_LR5_54HTB_440M'

ns={'bldg':'http://www.opengis.net/citygml/building/2.0','gml':'http://www.opengis.net/gml','gen':'http://www.opengis.net/citygml/generics/2.0'}


def pos3(text):
    vals=list(map(float,text.split())); pts=[tuple(vals[i:i+3]) for i in range(0,len(vals),3)]
    if pts and pts[0]==pts[-1]: pts=pts[:-1]
    return pts


def area_normal(pts):
    a=np.asarray(pts,float); n=np.zeros(3)
    for i in range(len(a)):
        x1,y1,z1=a[i]; x2,y2,z2=a[(i+1)%len(a)]
        n[0]+=(y1-y2)*(z1+z2); n[1]+=(z1-z2)*(x1+x2); n[2]+=(x1-x2)*(y1+y2)
    return 0.5*np.linalg.norm(n),n


def plane_basis(pts):
    a=np.asarray(pts,float); p0=a[0]; _,n=area_normal(pts)
    if np.linalg.norm(n)==0: return None
    if n[2]<0:n=-n
    n=n/np.linalg.norm(n)
    # longest plane-projected edge
    best=None; bl=0
    for i in range(len(a)):
        e=a[(i+1)%len(a)]-a[i]; e=e-np.dot(e,n)*n; L=np.linalg.norm(e)
        if L>bl: best=e; bl=L
    if best is None or bl<1e-9:
        t=np.array([1.,0,0]) if abs(n[0])<.9 else np.array([0,1.,0]); best=np.cross(n,t); bl=np.linalg.norm(best)
    u=best/bl; v=np.cross(n,u); v=v/np.linalg.norm(v)
    return p0,u,v,n


def xy_to_local(x,y,basis):
    p0,u,v,n=basis
    if abs(n[2])<1e-9: return None
    z=p0[2]-(n[0]*(x-p0[0])+n[1]*(y-p0[1]))/n[2]
    p=np.array([x,y,z]); d=p-p0
    return float(np.dot(d,u)),float(np.dot(d,v))


def plan_geom_to_local(g,basis):
    def one(poly):
        ext=[xy_to_local(x,y,basis) for x,y in list(poly.exterior.coords)]
        ext=[p for p in ext if p is not None]
        holes=[]
        for ring in poly.interiors:
            rr=[xy_to_local(x,y,basis) for x,y in list(ring.coords)]; rr=[p for p in rr if p is not None]
            if len(rr)>=4: holes.append(rr)
        return Polygon(ext, holes).buffer(0) if len(ext)>=4 else Polygon()
    if g.is_empty:return Polygon()
    if g.geom_type=='Polygon':return one(g)
    parts=[one(p) for p in g.geoms if p.geom_type=='Polygon']
    return unary_union([p for p in parts if not p.is_empty]) if parts else Polygon()


def dominant_angle(poly):
    mrr=poly.minimum_rotated_rectangle; c=list(mrr.exterior.coords); es=[]
    for a,b in zip(c[:-1],c[1:]):
        dx=b[0]-a[0];dy=b[1]-a[1]; es.append((math.hypot(dx,dy),math.degrees(math.atan2(dy,dx))))
    return max(es,key=lambda t:t[0])[1]


def _count_at(poly, angle, mw, ml, ox, oy):
    rp=rotate(poly,-angle,origin=poly.centroid,use_radians=False)
    minx,miny,maxx,maxy=rp.bounds; sx=mw+GAP; sy=ml+GAP
    xs=np.arange(minx+ox, maxx-mw+1e-9, sx)
    ys=np.arange(miny+oy, maxy-ml+1e-9, sy)
    if xs.size==0 or ys.size==0: return 0
    xx,yy=np.meshgrid(xs,ys)
    rects=shapely.box(xx.ravel(),yy.ravel(),xx.ravel()+mw,yy.ravel()+ml)
    return int(np.count_nonzero(shapely.covers(rp,rects)))


def pack(poly, mod, mode='standard'):
    if poly.is_empty:return 0
    parts=[poly] if poly.geom_type=='Polygon' else list(poly.geoms)
    total=0
    for p in parts:
        if p.area<mod['L']*mod['W']: continue
        a0=dominant_angle(p)
        if mode=='standard':
            angles=[a0,a0+90]; steps=3
        else:
            angles=list(np.arange(a0-15,a0+16,3))+list(np.arange(a0+75,a0+106,3)); steps=5
        best=0
        for ang in angles:
            # orientation is encoded by 90-degree angle set, module dimensions fixed
            mw,ml=mod['W'],mod['L']; sx=mw+GAP; sy=ml+GAP
            for ox in np.linspace(0,sx*.95,steps):
                for oy in np.linspace(0,sy*.95,steps):
                    best=max(best,_count_at(p,ang,mw,ml,ox,oy))
        total+=best
    return total


def parse_citygml():
    with zipfile.ZipFile(CITYGML_ZIP) as zf: root=etree.parse(zf.open('CityGML_Alles.gml')).getroot()
    rows=[]
    for b in root.xpath('.//bldg:Building',namespaces=ns):
        bid=b.get('{http://www.opengis.net/gml}id'); bname=''.join(b.xpath('./gml:name/text()',namespaces=ns)[:1])
        for p in b.xpath('./bldg:consistsOfBuildingPart/bldg:BuildingPart',namespaces=ns):
            pid=p.get('{http://www.opengis.net/gml}id'); attrs={a.get('name'): ''.join(a.xpath('./gen:value/text()',namespaces=ns)[:1]) for a in p.xpath('./gen:stringAttribute',namespaces=ns)}
            for rs in p.xpath('./bldg:boundedBy/bldg:RoofSurface',namespaces=ns):
                rid=rs.get('{http://www.opengis.net/gml}id')
                for j,pl in enumerate(rs.xpath('.//gml:Polygon/gml:exterior//gml:posList/text()',namespaces=ns),1):
                    pts=pos3(pl); plan0=Polygon([(x,y) for x,y,z in pts]) if len(pts)>=3 else Polygon()
                    was_valid=plan0.is_valid
                    plan=plan0 if was_valid else plan0.buffer(0)
                    area3d,n=area_normal(pts) if len(pts)>=3 else (0,np.zeros(3)); basis=plane_basis(pts) if len(pts)>=3 else None
                    key=hashlib.md5(plan.normalize().wkb + f"{min((q[2] for q in pts),default=0):.3f}|{max((q[2] for q in pts),default=0):.3f}".encode()).hexdigest() if not plan.is_empty else ''
                    rows.append({'building_id':bid,'building_name':bname,'part_id':pid,'roof_id':f'{rid}_{j}','alkisId':attrs.get('alkisId'),
                                 'npts':len(pts),'valid':was_valid,'repaired_valid':plan.is_valid,'area_plan':plan.area,'area_3d':area3d,'plan_wkt':plan.wkt,'pts3':pts,'basis':basis,'dupkey':key})
    return rows


def parse_alkis():
    ans=[]; tag='{http://www.adv-online.de/namespaces/adv/gid/6.0}AX_Gebaeude'; gml='{http://www.opengis.net/gml/3.2}'
    def p2(text):
        v=list(map(float,text.split())); return [tuple(v[i:i+2]) for i in range(0,len(v),2)]
    for _,e in etree.iterparse(str(ALKIS_XML),events=('end',),tag=tag,huge_tree=True):
        gid=e.get(gml+'id'); ring=e.find('.//'+gml+'Ring'); pts=[]
        if ring is not None:
            for pl in ring.xpath('.//gml:curveMember//gml:posList/text()',namespaces={'gml':'http://www.opengis.net/gml/3.2'}):
                seg=p2(pl)
                if not pts:pts.extend(seg)
                else:pts.extend(seg[1:] if pts[-1]==seg[0] else seg)
        if pts:
            if pts[0]!=pts[-1]:pts.append(pts[0])
            geom=Polygon(pts).buffer(0)
            if not geom.is_empty: ans.append((gid,geom))
        e.clear()
    return ans


def usable_local_for_row(row, building_envelope, setback):
    plan=wkt.loads(row['plan_wkt'])
    env=building_envelope.buffer(-setback, join_style=2) if setback>0 else building_envelope
    clipped=plan.intersection(env)
    return plan_geom_to_local(clipped,row['basis'])


def empirical_scenario(rows, setback, module_key, enhanced=False):
    mod=MODULES[module_key]; byb=defaultdict(list)
    for r in rows:byb[r['building_id']].append(r)
    facet_records=[]; building_records=[]
    for bid,rs in byb.items():
        plans=[wkt.loads(r['plan_wkt']) for r in rs]; envelope=unary_union(plans).buffer(0)
        plan_area=sum(p.area for p in plans); union_area=envelope.area; overlap=max(0,plan_area-union_area)
        # shared boundary length, exact boundaries
        shared=0.0
        for i in range(len(plans)):
            for j in range(i+1,len(plans)):
                shared += plans[i].boundary.intersection(plans[j].boundary).length
        areas=[]; comps=[]; usable_total=0; modules_total=0
        for r,p in zip(rs,plans):
            u=usable_local_for_row(r,envelope,setback); ua=u.area if not u.is_empty else 0.0
            n=pack(u,mod,'enhanced' if enhanced else 'standard') if ua>0 else 0
            areas.append(r['area_3d']); comps.append(4*np.pi*p.area/(p.length**2) if p.length>0 else np.nan)
            usable_total += ua; modules_total += n
            facet_records.append({'building_id':bid,'roof_id':r['roof_id'],'usable_area_m2':ua,'module_count':n,'roof_area_3d_m2':r['area_3d']})
        area_kwp=usable_total*(mod['kWp']/(mod['L']*mod['W'])); geom_kwp=modules_total*mod['kWp']
        ar=np.array(areas,float); small_thr=2*mod['L']*mod['W']; small_share=float(ar[ar<small_thr].sum()/ar.sum()) if ar.sum()>0 else np.nan
        building_records.append({'building_id':bid,'building_name':rs[0]['building_name'],'n_facets':len(rs),'parts':len(set(r['part_id'] for r in rs)),
            'sum_plan_area_m2':plan_area,'union_plan_area_m2':union_area,'plan_overlap_m2':overlap,'overlap_ratio':overlap/union_area if union_area else np.nan,
            'shared_boundary_m':shared,'internal_boundary_density_m_per_m2':shared/union_area if union_area else np.nan,
            'mean_facet_area_m2':float(ar.mean()),'median_facet_area_m2':float(np.median(ar)),'facet_area_cv':float(ar.std()/ar.mean()) if ar.mean()>0 else np.nan,
            'small_facet_area_share':small_share,'mean_compactness':float(np.nanmean(comps)),'usable_area_m2':usable_total,
            'facets_per_100m2':100*len(rs)/usable_total if usable_total else np.nan,'modules':modules_total,'geom_kwp':geom_kwp,'continuous_area_kwp':area_kwp,
            'packing_shortfall_pct':100*(1-geom_kwp/area_kwp) if area_kwp>0 else np.nan,
            'setback_m':setback,'module':module_key})
    return pd.DataFrame(building_records),pd.DataFrame(facet_records)


def alignment(rows, alkis):
    byb=defaultdict(list)
    for r in rows:byb[r['building_id']].append(wkt.loads(r['plan_wkt']))
    rec=[]
    for bid,plans in byb.items():
        env=unary_union(plans).buffer(0); c=env.centroid
        # candidate ALKIS intersecting env buffer 3m
        candidates=[(gid,g) for gid,g in alkis if g.intersects(env.buffer(3.0))]
        if candidates:
            au=unary_union([g for _,g in candidates])
            inter=env.intersection(au).area; union=env.union(au).area; iou=inter/union if union else np.nan
            # nearest centroid among candidate components
            cd=min(c.distance(g.centroid) for _,g in candidates)
            area_ratio=env.area/au.area if au.area else np.nan
            rec.append({'building_id':bid,'n_alkis_matches':len(candidates),'iou_roof_envelope_vs_alkis_union':iou,'centroid_offset_m':cd,'area_ratio_roof_to_alkis':area_ratio})
        else:
            rec.append({'building_id':bid,'n_alkis_matches':0,'iou_roof_envelope_vs_alkis_union':np.nan,'centroid_offset_m':np.nan,'area_ratio_roof_to_alkis':np.nan})
    return pd.DataFrame(rec)


def main():
    raw=parse_citygml(); rawdf=pd.DataFrame([{k:v for k,v in r.items() if k not in ('pts3','basis')} for r in raw])
    rawdf.to_csv(OUT/'citygml_raw_surface_audit.csv',index=False)
    # filter: valid, plan area > 0.05, no exact duplicates (none expected)
    keep=[]; seen=set(); exclusion=[]
    for r in raw:
        reason=None
        if (not r.get('repaired_valid',r['valid'])) or r['area_plan']<=0.05: reason='unrepairable_or_tiny_plan_polygon'
        elif (r['building_id'],r['dupkey']) in seen: reason='exact_duplicate'
        if reason: exclusion.append({'building_id':r['building_id'],'roof_id':r['roof_id'],'reason':reason,'area_plan':r['area_plan']})
        else: keep.append(r); seen.add((r['building_id'],r['dupkey']))
    pd.DataFrame(exclusion).to_csv(OUT/'citygml_exclusions.csv',index=False)
    print('raw surfaces',len(raw),'retained',len(keep),'excluded',len(exclusion),'buildings',len(set(r['building_id'] for r in keep)))

    # Base and sensitivity scenarios, external-envelope setback only.
    scenarios=[(0.0,BASE_MODULE),(0.3,BASE_MODULE),(0.5,BASE_MODULE),(0.3,'LONGi_LR5_72HGD_575M')]
    allsum=[]; base=None
    for sb,m in scenarios:
        print('scenario',sb,m,flush=True)
        b,f=empirical_scenario(keep,sb,m,enhanced=False)
        b.to_csv(OUT/f'empirical_buildings_{m}_setback_{sb:.1f}.csv',index=False)
        if sb==0.3 and m==BASE_MODULE: base=b; f.to_csv(OUT/'empirical_facets_base.csv',index=False)
        total_area=b.usable_area_m2.sum(); area_kwp=b.continuous_area_kwp.sum(); geom=b.geom_kwp.sum()
        allsum.append({'setback_m':sb,'module':m,'n_buildings':len(b),'usable_area_m2':total_area,'continuous_area_kwp':area_kwp,'packed_kwp':geom,
                       'aggregate_shortfall_pct':100*(1-geom/area_kwp),'median_building_shortfall_pct':b.packing_shortfall_pct.median(),
                       'q25':b.packing_shortfall_pct.quantile(.25),'q75':b.packing_shortfall_pct.quantile(.75)})
    pd.DataFrame(allsum).to_csv(OUT/'empirical_sensitivity_summary.csv',index=False)

    # fragmentation associations in base scenario
    vars=['n_facets','facets_per_100m2','internal_boundary_density_m_per_m2','median_facet_area_m2','small_facet_area_share','facet_area_cv','mean_compactness','usable_area_m2']
    corr=[]
    for v in vars:
        rho,p=spearmanr(base[v],base.packing_shortfall_pct,nan_policy='omit'); corr.append({'predictor':v,'spearman_rho':rho,'p_value':p})
    pd.DataFrame(corr).to_csv(OUT/'empirical_spearman.csv',index=False)

    # OLS with normalized fragmentation and size: boundary density + log usable area (standardized predictors)
    reg=base[['packing_shortfall_pct','internal_boundary_density_m_per_m2','usable_area_m2']].dropna().copy()
    reg['log_usable_area']=np.log(reg.usable_area_m2)
    for c in ['internal_boundary_density_m_per_m2','log_usable_area']:
        reg[c+'_z']=(reg[c]-reg[c].mean())/reg[c].std(ddof=0)
    X=sm.add_constant(reg[['internal_boundary_density_m_per_m2_z','log_usable_area_z']]); model=sm.OLS(reg.packing_shortfall_pct,X).fit(cov_type='HC3')
    pd.DataFrame({'term':model.params.index,'coef':model.params.values,'se_HC3':model.bse.values,'p_value':model.pvalues.values}).to_csv(OUT/'empirical_ols_HC3.csv',index=False)
    (OUT/'empirical_ols_summary.txt').write_text(model.summary().as_text())

    # leave-one-out largest areas
    loo=[]
    order=base.sort_values('usable_area_m2',ascending=False)
    for k in [0,1,3]:
        d=order.iloc[k:]
        loo.append({'removed_largest_n':k,'n_buildings':len(d),'aggregate_shortfall_pct':100*(1-d.geom_kwp.sum()/d.continuous_area_kwp.sum()),'packed_kwp':d.geom_kwp.sum(),'continuous_kwp':d.continuous_area_kwp.sum()})
    pd.DataFrame(loo).to_csv(OUT/'empirical_leave_largest_out.csv',index=False)

    # alignment audit
    alkis=parse_alkis(); al=alignment(keep,alkis); al.to_csv(OUT/'alkis_citygml_alignment.csv',index=False)
    print('alignment median IoU',al.iou_roof_envelope_vs_alkis_union.median(),'median centroid',al.centroid_offset_m.median(),'matched',sum(al.n_alkis_matches>0))

    # Algorithm benchmark on 12 moderate-sized facets across area quantiles using enhanced angle/offset search.
    # Reconstruct base usable local polygons and select facets in 10-250m2.
    byb=defaultdict(list)
    for r in keep: byb[r['building_id']].append(r)
    candidates=[]
    for bid,rs in byb.items():
        env=unary_union([wkt.loads(r['plan_wkt']) for r in rs]).buffer(0)
        for r in rs:
            u=usable_local_for_row(r,env,0.3)
            if 10<=u.area<=250:
                candidates.append((u.area,bid,r['roof_id'],u))
    candidates=sorted(candidates,key=lambda t:t[0])
    idx=np.linspace(0,len(candidates)-1,12).round().astype(int)
    bench=[]; mod=MODULES[BASE_MODULE]
    for ii in idx:
        a,bid,rid,u=candidates[ii]; std=pack(u,mod,'standard'); enh=pack(u,mod,'enhanced')
        bench.append({'building_id':bid,'roof_id':rid,'usable_area_m2':a,'standard_count':std,'enhanced_count':enh,'relative_gap_pct':100*(enh-std)/enh if enh>0 else 0})
        print('bench',rid,a,std,enh,flush=True)
    pd.DataFrame(bench).to_csv(OUT/'packing_algorithm_benchmark.csv',index=False)

    # compact summary
    summ={'raw_roof_polygons':len(raw),'retained_roof_polygons':len(keep),'excluded':len(exclusion),'buildings':base.shape[0],
          'base_usable_area_m2':float(base.usable_area_m2.sum()),'base_continuous_kwp':float(base.continuous_area_kwp.sum()),'base_packed_kwp':float(base.geom_kwp.sum()),
          'base_aggregate_shortfall_pct':float(100*(1-base.geom_kwp.sum()/base.continuous_area_kwp.sum())),
          'base_median_shortfall_pct':float(base.packing_shortfall_pct.median()),'base_iqr':[float(base.packing_shortfall_pct.quantile(.25)),float(base.packing_shortfall_pct.quantile(.75))],
          'alignment_median_iou':float(al.iou_roof_envelope_vs_alkis_union.median()),'alignment_median_centroid_offset_m':float(al.centroid_offset_m.median()),
          'algorithm_benchmark_median_gap_pct':float(pd.DataFrame(bench).relative_gap_pct.median()),'algorithm_benchmark_max_gap_pct':float(pd.DataFrame(bench).relative_gap_pct.max())}
    (OUT/'hardening_summary.json').write_text(json.dumps(summ,indent=2))
    print(json.dumps(summ,indent=2))

if __name__=='__main__':
    raise SystemExit('Internal implementation module. Use run_scenario.py, run_benchmark.py, or reproduce_empirical.py.')
