"""Compile the manuscript with Tectonic into an explicit external build directory.

The compiler is supplied by the caller. Compilation products stay in the
selected build directory; final PDF and diagnostic records are copied into
this repository only for the default paper source.
"""
from pathlib import Path
import argparse
import calendar
import hashlib
import json
import os
import re
import shutil
import subprocess
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parents[2]
AUDIT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--source', type=Path, default=ROOT/'paper')
    parser.add_argument('--cache', type=Path)
    parser.add_argument('--archive-check', action='store_true')
    args = parser.parse_args()
    source = args.source.resolve()
    out = args.output_dir.resolve()
    assert not out.is_relative_to(source), 'Build outside the manuscript source'
    assert args.compiler.is_file() and (source/'main.tex').is_file()
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env['SOURCE_DATE_EPOCH'] = str(calendar.timegm((2026,11,3,0,0,0)))
    if args.cache:
        env['TECTONIC_CACHE_DIR'] = str(args.cache.resolve())
    if os.name == 'nt' and 'WINDIR' in env:
        conf = out/'fonts.conf'
        fonts = (Path(env['WINDIR'])/'Fonts').as_posix()
        conf.write_text('<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "fonts.dtd">'
                        '<fontconfig><dir>'+escape(fonts)+'</dir><cachedir>'
                        +escape((out/'fontcache').as_posix())+'</cachedir></fontconfig>',
                        encoding='utf-8')
        env['FONTCONFIG_FILE'] = str(conf)
    command = [str(args.compiler.resolve()), '--keep-logs', '--keep-intermediates',
               '--outdir', str(out), str(source/'main.tex')]
    result = subprocess.run(command, cwd=source, env=env, text=True,
                            encoding='utf-8', errors='replace',
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    prefix = 'archive_' if args.archive_check else ''
    (AUDIT/(prefix+'compile.log')).write_text(result.stdout, encoding='utf-8')
    for suffix in ['log','bbl']:
        path = out/('main.'+suffix)
        if path.exists():
            shutil.copy2(path, AUDIT/(prefix+'main.'+suffix))
    log = (out/'main.log').read_text(encoding='utf-8', errors='replace') if (out/'main.log').exists() else ''
    diagnostics = {
        'revision': 'v1.0.1', 'returncode': result.returncode,
        'source_sha256': hashlib.sha256((source/'main.tex').read_bytes()).hexdigest(),
        'compiler': str(args.compiler.resolve()),
        'requested_warning_counts': {
            'microtype_footnote_patch': len(re.findall(r'Unable to apply patch.*footnote', log)),
            'underfull_hbox': len(re.findall(r'Underfull \\hbox', log)),
        },
        'overfull_box_count': len(re.findall(r'Overfull \\[hv]box', log)),
        'unresolved_reference_count': len(re.findall(r'(?:Citation|Reference).*undefined|There were undefined references',log)),
        'engine_warnings': [line for line in result.stdout.splitlines() if 'warning:' in line.lower()],
        'archive_check': args.archive_check,
    }
    if result.returncode:
        print(result.stdout[-6000:])
        raise SystemExit(result.returncode)
    assert (out/'main.pdf').is_file()
    if not args.archive_check:
        shutil.copy2(out/'main.pdf', ROOT/'manuscript.pdf')
    (AUDIT/(prefix+'compilation.json')).write_text(json.dumps(diagnostics,indent=2),encoding='utf-8')
    print(json.dumps(diagnostics,indent=2))


if __name__ == '__main__':
    main()
