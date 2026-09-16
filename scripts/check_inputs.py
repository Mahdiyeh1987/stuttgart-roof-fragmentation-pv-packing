from pathlib import Path
import hashlib,pandas as pd
root=Path(__file__).resolve().parents[1]
expected=pd.read_csv(root/'results'/'input_file_manifest_sha256.csv')
fail=[]
for _,r in expected.iterrows():
 p=root/'data_raw'/r.filename
 h=hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
 ok=(h==r.sha256 and p.stat().st_size==r.bytes) if p.exists() else False
 print(r.filename,'OK' if ok else 'MISMATCH')
 if not ok: fail.append(r.filename)
raise SystemExit(1 if fail else 0)
