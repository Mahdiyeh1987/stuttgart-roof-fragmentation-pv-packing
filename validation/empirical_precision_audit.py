"""Controlled 10-nm numerical-tolerance audit, not a physically larger usable roof."""
from pathlib import Path
import sys,time,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import analysis_core as core
import numpy as np, pandas as pd, shapely
from shapely.affinity import rotate

def robust_count_angle_gap(poly,angle,mod,gap,steps):
 rp=rotate(poly,-angle,origin=poly.centroid,use_radians=False)
 pbuf=rp.buffer(1.e-8) # numerical precision: 10 nanometres, not rooftop setback
 minx,miny,maxx,maxy=rp.bounds
 mw,ml=mod['W'],mod['L'];sx=mw+gap;sy=ml+gap
 best=0
 for ox in np.linspace(0,sx*.98,steps):
  xs=np.arange(minx+ox,maxx-mw+1e-9,sx)
  if not xs.size:continue
  for oy in np.linspace(0,sy*.98,steps):
   ys=np.arange(miny+oy,maxy-ml+1e-9,sy)
   if not ys.size:continue
   xx,yy=np.meshgrid(xs,ys)
   rects=shapely.box(xx.ravel(),yy.ravel(),xx.ravel()+mw,yy.ravel()+ml)
   best=max(best,int(np.count_nonzero(shapely.covers(pbuf,rects))))
 return best

def run(outdir):
 outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
 t0=time.monotonic()
 core.count_angle_gap=robust_count_angle_gap
 raw,k0,_=core.initial_rows();k1,_,_=core.remove_nested(k0); facets,_,_=core.merge_coplanar_adjacent(k1)
 print('facets',len(facets),'prep seconds',round(time.monotonic()-t0,2),flush=True)
 b,f=core.scenario_fast(facets,setback=.3,mkey=core.BASE,gap=.01,workers=4)
 old=pd.read_csv(Path(__file__).resolve().parents[1]/'results'/'base_facets.csv')
 delta=f.merge(old[['physical_facet_uid','modules']],on='physical_facet_uid',suffixes=('_robust','_archive'))
 delta['delta_modules']=delta.modules_robust-delta.modules_archive
 print('n compare',len(delta),'differences',(delta.delta_modules!=0).sum(), 'total delta',delta.delta_modules.sum(),'min delta',delta.delta_modules.min(),'max delta',delta.delta_modules.max(),flush=True)
 orig=pd.read_csv(Path(__file__).resolve().parents[1]/'results'/'base_buildings.csv')
 orig_short=100*(1-orig.packed_kwp.sum()/orig.continuous_kwp.sum())
 new_short=100*(1-b.packed_kwp.sum()/b.continuous_kwp.sum())
 print('original',orig_short,'robust',new_short,'elapsed',round(time.monotonic()-t0,2),flush=True)
 b.to_csv(outdir/'precision_buffer_buildings.csv',index=False);f.to_csv(outdir/'precision_buffer_facets.csv',index=False)
 delta.to_csv(outdir/'precision_buffer_facet_comparison.csv',index=False)
 s={'original_shortfall_pct':orig_short,'precision_tolerant_shortfall_pct':new_short,'delta_shortfall_pp':new_short-orig_short,'n_physical_facets':len(f),'facets_changed':int((delta.delta_modules!=0).sum()),'module_count_difference_total':int(delta.delta_modules.sum()),'runtime_s':round(time.monotonic()-t0,2),'numerical_buffer_m':1e-8}
 (outdir/'precision_audit_summary.json').write_text(json.dumps(s,indent=2));return s
if __name__=='__main__':print(run(Path(__file__).parent/'outputs'),flush=True)
