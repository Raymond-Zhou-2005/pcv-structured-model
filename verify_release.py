"""Verify the downloaded release before regenerating its scientific outputs."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def included(path):
    relative = path.relative_to(ROOT)
    return (path.is_file() and path.name != 'release_manifest.json'
            and not any(part in {'.git', '__pycache__', '.venv'} for part in relative.parts)
            and not relative.as_posix().startswith('checks/logs/')
            and path.suffix not in {'.zip', '.pyc', '.pyo'})


def main():
    manifest = json.loads((ROOT / 'release_manifest.json').read_text(encoding='utf-8'))
    expected = {row['path']: row for row in manifest['files']}
    actual = {p.relative_to(ROOT).as_posix(): p for p in ROOT.rglob('*') if included(p)}
    assert set(expected) == set(actual), {
        'missing': sorted(set(expected) - set(actual)),
        'unexpected': sorted(set(actual) - set(expected))}
    for relative, record in expected.items():
        payload = actual[relative].read_bytes()
        assert len(payload) == record['bytes'], relative
        assert hashlib.sha256(payload).hexdigest() == record['sha256'], relative
    for relative in ['LICENSE', 'LICENSE-DATA', 'LICENSES.md', 'CITATION.cff',
                     'docs/artifact_map.csv', 'docs/data_dictionary.md',
                     'docs/availability.md', 'manuscript.pdf', 'paper/main.tex']:
        assert relative in expected, relative
    print(json.dumps({'release': manifest['release'], 'verified_files': len(expected),
                      'status': 'passed',
                      'note': 'Run before reproduction; generated outputs can change within declared tolerances.'}, indent=2))


if __name__ == '__main__':
    main()
