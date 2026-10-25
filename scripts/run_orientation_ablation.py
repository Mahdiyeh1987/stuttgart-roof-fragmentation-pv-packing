"""2026 revision-only fixed-axis orientation ablation and runtime check.

This is a NEW analysis, separate from the archived 2026 original study.
Compares module grids oriented at the dominant roof-plane axis vs 90deg offset,
using identical geometry, external setback, module and 9x9 search.
Does not represent rack tilt or true flat-roof deployability.
"""
from __future__ import annotations
import os, math, time, sys, platform
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
from shapely import wkt, set_precision
from shapely.ops import unary_union
import analysis_core as v
import hardening_analysis as ha

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'revised_checks'
OUT.mkdir(exist_ok=True)

def _task(task):
    bid, fid, poly_wkt, orientation, angle, modkey, gap=task
    poly=wkt.loads(poly_wkt)
    mod=ha.MODULES[modkey]
    n=v.count_angle_gap(poly,angle,mod,gap,9)
    return {'building_id':bid, 'physical_facet_uid':fid, 'part_orientation':orientation, 'modules':int(n)}

def main():
    t0=time.perf_counter()
    raw,ret,_=v.initial_rows()
    kept,_,_=v.remove_nested(ret)
    facets,_,_=v.merge_coplanar_adjacent(kept)
    prepare_s=time.perf_counter()-t0
    by=defaultdict(list)
    for facet in facets:by[facet['building_id']].append(facet)
    task=[]; geomrows=[]
    m=ha.MODULES[v.BASE]
    for bid,fs in by.items():
        plans=[set_precision(wkt.loads(f['plan_wkt']),v.GRID) for f in fs]
        env=wkt.loads(unary_union(plans).buffer(0).wkt)
        e=env.buffer(-.30, join_style=2)
        for f in fs:
            clipped=wkt.loads(f['plan_wkt']).intersection(e)
            u=ha.plan_geom_to_local(clipped,f['basis'])
            parts=([u] if u.geom_type=='Polygon' else [p for p in getattr(u,'geoms',[]) if p.geom_type=='Polygon']) if not u.is_empty else []
            for i,p in enumerate(parts):
                key=f"{f['physical_facet_uid']}__{i}"
                geomrows.append({'building_id':bid,'physical_facet_uid':key,'usable_area_m2':float(p.area)})
                if p.area<m['L']*m['W']:continue
                a=ha.dominant_angle(p)%180
                for label,ang in [('axis_0',a),('axis_90',(a+90)%180)]:
                    task.append((bid,key,p.wkt,label,float(ang),v.BASE,.01))
    candidates=[]
    packing_t0=time.perf_counter()
    with ProcessPoolExecutor(max_workers=min(5,os.cpu_count() or 1)) as ex:
        fs=[ex.submit(_task,t) for t in task]
        for f in as_completed(fs):candidates.append(f.result())
    packing_s=time.perf_counter()-packing_t0
    cd=pd.DataFrame(candidates); geometries=pd.DataFrame(geomrows)
    pv=cd.pivot_table(index=['building_id','physical_facet_uid'],columns='part_orientation',values='modules',aggfunc='max').fillna(0).reset_index()
    d=geometries.merge(pv,on=['building_id','physical_facet_uid'],how='left').fillna(0)
    d['axis_best']=d[['axis_0','axis_90']].max(axis=1)
    d['continuous_kwp']=d.usable_area_m2*m['kWp']/(m['L']*m['W'])
    for mode in ['axis_0','axis_90','axis_best']:
        modcount=d[mode].sum(); cont=d.continuous_kwp.sum(); packed=modcount*m['kWp'];short=100*(1-packed/cont)
        print(f'{mode} modules={int(modcount)} capacity_kwp={packed:.3f} continuous_kwp={cont:.3f} shortfall_pct={short:.7f}',flush=True)
    print(f'raw={len(raw)} retained={len(kept)} facets={len(facets)} candidate_tasks={len(task)}')
    print(f'prep_sec={prepare_s:.3f} packing_sec={packing_s:.3f} elapsed_sec={time.perf_counter()-t0:.3f}')
    print(f'python={platform.python_version()} cpu_count={os.cpu_count()} method=ProcessPoolExecutor workers={min(5,os.cpu_count() or 1)}')
    d.to_csv(OUT/'orientation_axis_ablation_facets.csv',index=False)
    cd.to_csv(OUT/'orientation_axis_ablation_candidates.csv',index=False)
    msummary={'method':'two orthogonal dominant-axis placements with 9x9 offset search; not an unrestricted layout or rack-tilt study',
       'n_buildings':len(by),'n_physical_facets':len(facets),'n_parts':len(d),'workers':min(5,os.cpu_count() or 1),
       'preparation_seconds':prepare_s,'packing_seconds':packing_s,'total_seconds':time.perf_counter()-t0,'python_version':platform.python_version(),
       'continuous_kwp':float(d.continuous_kwp.sum()),
       'orientation_modes':{mode:{'n_modules':int(d[mode].sum()),'packed_kwp':float(d[mode].sum()*m['kWp']), 'packing_shortfall_pct':float(100*(1-d[mode].sum()*m['kWp']/d.continuous_kwp.sum()))} for mode in ['axis_0','axis_90','axis_best']}}
    import json
    (OUT/'orientation_axis_ablation_summary.json').write_text(json.dumps(msummary,indent=2))
if __name__=='__main__':main()
