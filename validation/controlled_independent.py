"""Independent (not historical) controlled-subdivision replication.
Uses fixed roof exterior after 0.30 m setback and same archived deterministic packing routine.
This DOES NOT recover historical random polygons and does not guarantee numerical identity.
"""
from pathlib import Path
import sys, argparse, json, time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import analysis_core as a
import numpy as np, pandas as pd
from shapely.geometry import box,LineString
from shapely.ops import split, unary_union
from scipy.stats import spearmanr

SEED=20260913
MIN_PIECE_M2=.8
MOD=a.ha.MODULES[a.BASE]
ROOT=box(.30,.30,11.70,9.70)

def random_point_inside(poly,rng):
 from shapely.geometry import Point
 x0,y0,x1,y1=poly.bounds
 for _ in range(3000):
  pt=Point(rng.uniform(x0,x1),rng.uniform(y0,y1))
  if poly.contains(pt):return pt
 raise RuntimeError('point-rejection sampler exceeded limit')

def split_one_piece(p,rng):
 point=random_point_inside(p,rng)
 theta=rng.uniform(0,np.pi)
 direction=np.array([np.cos(theta),np.sin(theta)])
 D=2000
 ln=LineString([(point.x-D*direction[0],point.y-D*direction[1]),(point.x+D*direction[0],point.y+D*direction[1])])
 pieces=[g for g in split(p,ln).geoms if g.area>1e-10]
 if len(pieces)!=2 or min(g.area for g in pieces)<MIN_PIECE_M2:return None
 if abs(sum(g.area for g in pieces)-p.area)>1e-6: raise RuntimeError('area conservation failed')
 return pieces

def construct_fragments(n,rng):
 parts=[ROOT]
 for _ in range(n-1):
  for trial in range(5000):
   probs=np.array([p.area for p in parts]); probs/=probs.sum()
   ix=int(rng.choice(len(parts),p=probs))
   pieces=split_one_piece(parts[ix],rng)
   if pieces is not None:
    parts=parts[:ix]+pieces+parts[ix+1:]
    break
  else:raise RuntimeError(f'No valid split for n={n}')
 assert len(parts)==n
 assert abs(sum(p.area for p in parts)-107.16)<1e-5
 assert len([p for p in parts if p.area<MIN_PIECE_M2-1e-9])==0
 return parts

def run(repeats=50,seed=SEED,outdir=None):
 rng=np.random.default_rng(seed)
 n0=a.pack_production(ROOT,MOD,.01)
 print('base_count=',n0,'area=',ROOT.area,flush=True)
 if n0!=48: raise RuntimeError('unfragmented N0 differs from article 48')
 rec=[];t0=time.monotonic()
 for n in range(2,13):
  for rep in range(repeats):
   parts=construct_fragments(n,rng)
   counts=[a.pack_production(p,MOD,.01) for p in parts]
   area=np.array([p.area for p in parts]); boundary=sum(p.length for p in parts)-ROOT.length
   compact=np.mean([4*np.pi*p.area/(p.length**2) for p in parts])
   small=float(area[area<2*MOD['L']*MOD['W']].sum()/area.sum())
   rec.append(dict(facet_count=n,replicate=rep+1,n_modules=sum(counts),incremental_loss_pct=100*(1-sum(counts)/n0),total_area_m2=area.sum(),internal_boundary_density=boundary/(2*ROOT.area),small_facet_area_share=small,facet_area_cv=area.std(ddof=0)/area.mean(),mean_compactness=compact))
  print('facets',n,'median loss',float(np.median([x['incremental_loss_pct'] for x in rec if x['facet_count']==n])),'runtime_sec',round(time.monotonic()-t0,1),flush=True)
 df=pd.DataFrame(rec)
 summary=df.groupby('facet_count').incremental_loss_pct.agg(['median',lambda x:np.quantile(x,.25),lambda x:np.quantile(x,.75),'min','max']).reset_index();summary.columns=['facet_count','median_loss_pct','p25_pct','p75_pct','min_pct','max_pct']
 historical=pd.read_csv(Path(__file__).resolve().parents[1]/'results'/'controlled_production_summary.csv')
 print('HISTORICAL cols',historical.columns.tolist());print(historical.to_string(index=False)[:2000],flush=True)
 correl=[]
 for col in ['facet_count','internal_boundary_density','small_facet_area_share','facet_area_cv','mean_compactness']:
  stat=spearmanr(df[col],df.incremental_loss_pct)
  correl.append(dict(variable=col,spearman_rho=stat.statistic,p=stat.pvalue))
 print(summary.to_string(index=False),flush=True);print(pd.DataFrame(correl).to_string(index=False),flush=True)
 if outdir:
  outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
  df.to_csv(outdir/'independent_controlled_realisations.csv',index=False)
  summary.to_csv(outdir/'independent_controlled_summary.csv',index=False)
  pd.DataFrame(correl).to_csv(outdir/'independent_controlled_correlations.csv',index=False)
  (outdir/'metadata.json').write_text(json.dumps(dict(seed=seed,repeats=repeats,outer_roof_m=[12,10],external_setback_m=.3,min_piece_area_m2=.8,uniform_interior_point_sampler='rejection bounding box',split_angle='uniform [0,pi)',area_weighted_facet_choice=True,rejection_policy='retry splits creating part below 0.8 m²',same_pack_production=True,independent_generator_not_identical_to_historical=True,baseline_modules=n0,seconds_elapsed=round(time.monotonic()-t0,3)),indent=2))
 return df,summary

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--repeats',type=int,default=50);parser.add_argument('--seed',type=int,default=SEED);parser.add_argument('--outdir',default=str(Path(__file__).parent/'outputs'))
 args=parser.parse_args();run(args.repeats,args.seed,args.outdir)
