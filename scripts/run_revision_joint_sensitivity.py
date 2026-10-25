"""Additional 3x2 setback and gap factorial sensitivity; 2026 revision.
Uses exact original empirical algorithm, not a simulation of rooftop installation.
"""
import os,sys,json,time,platform
from pathlib import Path
import pandas as pd
import analysis_core as v

def main():
 root=Path(__file__).resolve().parents[1]
 out=root/'revised_checks';out.mkdir(exist_ok=True)
 t=time.perf_counter()
 raw,ret,_=v.initial_rows()
 kept,_,_=v.remove_nested(ret)
 f,_,_=v.merge_coplanar_adjacent(kept)
 rows=[]
 for sb in (float(sys.argv[1]),):
  for gap in (float(sys.argv[2]),):
   t0=time.perf_counter()
   b,ff=v.scenario_fast(f, sb, v.BASE, gap, None, workers=min(5,os.cpu_count() or 1))
   p=float(b.packed_kwp.sum());r=float(b.continuous_kwp.sum());s=100*(1-p/r)
   rec={'setback_m':sb,'clearance_mm':int(gap*1000),'n_buildings':len(b),'n_modules':int(b.modules.sum()),'usable_area_m2':float(b.usable_area_m2.sum()),'continuous_kwp':r,'packed_kwp':p,'shortfall_pct':s,'runtime_seconds':time.perf_counter()-t0}
   rows.append(rec)
   print(json.dumps(rec),flush=True)
 csv=out/'setback_clearance_joint_sensitivity.csv'
 old=pd.read_csv(csv) if csv.exists() else pd.DataFrame()
 pd.concat([old,pd.DataFrame(rows)],ignore_index=True).drop_duplicates(['setback_m','clearance_mm'],keep='last').to_csv(csv,index=False)
 info={'preparation_plus_scenarios_total_seconds':time.perf_counter()-t,'python':platform.python_version(),'cpu_count':os.cpu_count(),'method':'same production search, 9x9 offsets, all 340 facets and 31 buildings','n_input_facets':len(f)}
 (out/'setback_clearance_joint_sensitivity_metadata.json').write_text(json.dumps(info,indent=2))
 print('METADATA',info,flush=True)
if __name__=='__main__':main()
