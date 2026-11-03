"""Write a non-circular, sorted SHA-256 manifest for the repository release."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent

def main():
    files=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(ROOT)
        if p.name=='release_manifest.json' or any(x in rel.parts for x in ['.git','__pycache__','.venv']):continue
        if rel.as_posix().startswith('checks/logs/') or p.suffix in ['.zip','.pyc','.pyo']:continue
        files.append(dict(path=rel.as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    result=dict(release='v1.0.1',repository='https://github.com/Raymond-Zhou-2005/pcv-structured-model',
                release_date='2026-11-03',publication_status='versioned public research release',
                base_release='v1.0.0',base_commit='6e9885d8c605a10fb2e73a27288c2d48d3968cdf',
                manifest_scope='All retained release files except this manifest, transient logs, caches and nested ZIPs',
                evidence_scope='Arithmetic and computational reproducibility; empirical parameter qualifications remain in the ledgers',
                files=files)
    (ROOT/'release_manifest.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('Manifest:',len(files),'files')

if __name__=='__main__':main()
