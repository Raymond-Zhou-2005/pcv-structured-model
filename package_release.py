"""Package a verified repository and the compiled manuscript for a GitHub release.

Compile paper/main.tex first; place its PDF at manuscript.pdf. Run the numerical
checks and write_manifest.py before this script. Output must be outside the
repository to prevent recursive/nested release files.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def archive(destination, members):
    assert not destination.exists(), destination
    with zipfile.ZipFile(destination, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for relative, path in sorted(members.items()):
            info = zipfile.ZipInfo(relative, (2026, 10, 28, 0, 0, 0))
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            z.writestr(info, path.read_bytes(), compress_type=zipfile.ZIP_DEFLATED,
                       compresslevel=6)
    with zipfile.ZipFile(destination) as z:
        assert z.testzip() is None
        for relative, path in members.items():
            assert hashlib.sha256(z.read(relative)).hexdigest() == digest(path), relative


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    assert not output.is_relative_to(ROOT), 'Release archives must be outside the repository'
    output.mkdir(parents=True, exist_ok=True)
    subprocess.run([sys.executable, str(ROOT / 'verify_release.py')], check=True)
    checks = json.loads((ROOT / 'checks/verification.json').read_text())
    assert checks['passed'] == checks['check_count'] and not checks['failed']
    paper = {p.relative_to(ROOT / 'paper').as_posix(): p
             for p in (ROOT / 'paper').rglob('*') if p.is_file()}
    assert len(paper) == 9, 'Only the nine actual Overleaf source dependencies are expected'
    archive(output / 'overleaf.zip', paper)
    manifest = json.loads((ROOT / 'release_manifest.json').read_text())
    payload = {row['path']: ROOT / row['path'] for row in manifest['files']}
    payload['release_manifest.json'] = ROOT / 'release_manifest.json'
    archive(output / 'github_repository.zip', payload)
    for name in ['manuscript.pdf', 'README.md']:
        destination = output / name
        assert not destination.exists(), destination
        destination.write_bytes((ROOT / name).read_bytes())
    archive(output / 'PCV_Main_Revision.zip', {
        name: output / name for name in
        ['README.md', 'manuscript.pdf', 'overleaf.zip', 'github_repository.zip']})
    assets = [dict(name=name, bytes=(output / name).stat().st_size,
                   sha256=digest(output / name)) for name in
              ['manuscript.pdf', 'overleaf.zip', 'github_repository.zip', 'PCV_Main_Revision.zip']]
    (output / 'SHA256SUMS.txt').write_text(
        ''.join(f"{row['sha256']}  {row['name']}\n" for row in assets), encoding='ascii')
    (output / 'release_assets.json').write_text(json.dumps(assets, indent=2), encoding='utf-8')
    print(json.dumps({'assets': assets, 'repository_files': len(payload)}, indent=2))


if __name__ == '__main__':
    main()
