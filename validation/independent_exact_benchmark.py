"""Independent analytic-oracle checks; no historical verification claim.
Oracle = geometric construction + area upper bound, except narrow strip proven by width (< module minimum).
Tests cover rotations, setbacks/margins, concavity, hole, clearance, width exclusion.
"""
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import analysis_core as core
import pandas as pd
from shapely.geometry import box
from shapely.ops import unary_union
from shapely.affinity import rotate
from empirical_precision_audit import robust_count_angle_gap

M=core.ha.MODULES[core.BASE];w,l=M['W'],M['L']
def correct(poly,gap):
 res=0
 for part in ([poly] if poly.geom_type=='Polygon' else list(poly.geoms)):
  if part.area+1e-10<w*l:continue
  res+=max(robust_count_angle_gap(part,float(angle),M,gap,9) for angle in core.edge_angles(part))
 return int(res)

def main(out_dir):
 out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
 data=[]
 def check(id,poly,expected,gap=0,proof='explicit tiling + area upper bound'):
  raw=core.pack_production(poly,M,gap)
  fix=correct(poly,gap)
  data.append(dict(case=id,expected_exact_modules=expected,production_modules=raw,tolerance_corrected_modules=fix,production_pass=(raw==expected),corrected_pass=(fix==expected),module_gap_m=gap,roof_area_m2=round(poly.area,9),area_upper_bound=int(poly.area/(w*l)+1e-8),proof=proof))
 for a,b in [(1,1),(2,2),(3,2),(4,3)]:check(f'axis_exact_{a}x{b}',box(0,0,a*w,b*l),a*b)
 for angle in [23,47.5,90]:check(f'rotated_exact_{angle}',rotate(box(0,0,3*w,2*l),angle,origin='centroid'),6)
 t=lambda x,y:box(x*w,y*l,(x+1)*w,(y+1)*l)
 check('L_shaped_three_modules',unary_union([t(0,0),t(1,0),t(0,1)]),3)
 check('roof_hole_five_modules',box(0,0,3*w,2*l).difference(t(1,0)),5)
 check('narrow_straight_strip',box(0,0,.7,8),0,proof='width 0.7 m < minimum 1.134 m dimension, for any rotation')
 for a,b in [(1,1),(2,2),(3,2)]:check(f'clearance_10mm_exact_{a}x{b}',box(0,0,a*w+(a-1)*.01,b*l+(b-1)*.01),a*b,.01,proof='constructive lattice + area bound')
 for angle in [23,47.5]:
  check(f'rotated_5cm_margin_{angle}',rotate(box(-.05,-.05,3*w+.05,2*l+.05),angle,origin='centroid'),6)
 for angle in [23,47.5]:
  check(f'rotated_clearance_10mm_margin_{angle}',rotate(box(-.05,-.05,3*w+2*.01+.05,2*l+.01+.05),angle,origin='centroid'),6,.01)
 df=pd.DataFrame(data);df.to_csv(out_dir/'exact_oracle_validation_cases.csv',index=False)
 summary=dict(test_cases=len(df),passed_original=int(df.production_pass.sum()),failed_original=int((~df.production_pass).sum()),passed_precision_corrected=int(df.corrected_pass.sum()),failed_precision_corrected=int((~df.corrected_pass).sum()),expected_values_are_proven_by_constructive_tiling_and_area_bound=True,numerical_tolerance_m=1e-8)
 (out_dir/'exact_oracle_summary.json').write_text(json.dumps(summary,indent=2))
 print(df[['case','expected_exact_modules','production_modules','tolerance_corrected_modules','production_pass','corrected_pass']].to_string(index=False));print(summary)
 return df
if __name__=='__main__':main(Path(__file__).parent/'outputs')
