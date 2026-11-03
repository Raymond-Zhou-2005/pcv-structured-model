"""Validate revised text, figures, compilation and unchanged numerical assets."""
from pathlib import Path
import argparse
import ast
import collections
import hashlib
import json
import re

import fitz
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
AUDIT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def citations(text):
    return collections.Counter(k.strip() for group in re.findall(r'\\cite\{([^}]+)\}',text)
                               for k in group.split(','))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    args = parser.parse_args()
    checks = []

    def check(name, passed, detail=''):
        checks.append(dict(check=name, passed=bool(passed), detail=detail))

    text = (ROOT/'paper/main.tex').read_text(encoding='utf-8')
    bib = (ROOT/'paper/references.bib').read_text(encoding='utf-8')
    keys = set(re.findall(r'@\w+\{([^,]+)',bib))
    used = citations(text)
    check('citation_keys_resolve',set(used).issubset(keys),sorted(set(used)-keys))
    check('original_24_citation_keys',len(used)==24,len(used))
    check('original_25_citation_occurrences',sum(used.values())==25,sum(used.values()))
    for marker in ['Author-review version','require author confirmation',
                   'Acknowledgment wording','CFR_TRANSLATION_AUDIT_RESULT']:
        check('removed_placeholder',marker not in text,marker)
    check('AI_disclosure_not_added',not re.search(r'AI use|ChatGPT|Codex|large language model',text,re.I))
    check('interval_density_counterexample',all(x in text for x in ['eq:interval_cross_difference','u_i=\\mathbf1_{E_i}/|E_i|']))
    check('local_positive_denominator_tube',all(x in text for x in ['D_0=D(A_0)>0','eq:local_contraction']))
    check('CFR_translation_formula_present','eq:cfr_translation' in text)
    check('CFR_frozen_exposure_quantification','0.014--0.081\\%' in text and 'exposure fixed' in text)
    check('microtype_footnote_patch_disabled','\\microtypesetup{nopatch=footnote}' in text)
    check('bibliography_local_alignment','\\AtBeginEnvironment{thebibliography}{\\raggedright}' in text)
    check('warnings_not_hidden',not re.search(r'\\(?:hbadness|hfuzz|sloppy)\b',text))
    compilation = json.loads((AUDIT/'compilation.json').read_text(encoding='utf-8'))
    check('compiled_current_source',compilation['source_sha256']==sha(ROOT/'paper/main.tex'))
    check('compiler_returncode',compilation['returncode']==0)
    for key,value in compilation['requested_warning_counts'].items():
        check('requested_warning_removed',value==0,key+': '+str(value))
    check('no_overfull_boxes',compilation['overfull_box_count']==0,compilation['overfull_box_count'])
    check('no_unresolved_references',compilation['unresolved_reference_count']==0)
    with fitz.open(ROOT/'manuscript.pdf') as doc:
        all_text='\n'.join(p.get_text() for p in doc)
        check('PDF_no_question_mark_crossrefs','??' not in all_text)
        check('PDF_title_present','From cohort state transitions' in all_text)
        check('PDF_all_six_figure_captions',all(re.search(r'Fig\s*'+str(i)+r'\s*:',all_text) for i in range(1,7)))
        check('PDF_no_replacement_characters','\ufffd' not in all_text)
        outside=[]
        for i,page in enumerate(doc):
            for block in page.get_text('dict')['blocks']:
                for line in block.get('lines',[]):
                    for span in line.get('spans',[]):
                        x0,y0,x1,y1=span['bbox']
                        if x0 < -.5 or y0 < -.5 or x1 > page.rect.width+.5 or y1 > page.rect.height+.5:
                            outside.append(dict(page=i+1,text=span['text'],bbox=span['bbox']))
        check('PDF_no_content_outside_page',not outside,outside)
        page_count=len(doc)
    for n in range(1,7):
        pdf=ROOT/'figures'/f'Fig{n}.pdf'
        check('paper_figure_matches_export',sha(pdf)==sha(ROOT/'paper/figures'/pdf.name),n)
        with fitz.open(pdf) as doc:
            page=doc[0]
            spans=[s for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l.get('spans',[]) if s['text'].strip()]
            check('figure_dimensions',page.rect.width<=7.5*72+.01 and page.rect.height<=8.75*72+.01,n)
            check('figure_fonts_8_to_12',all(7.99<=s['size']<=12.01 for s in spans),n)
            check('figure_allowed_font_families',all(any(k in s['font'] for k in ['Arial','Times','Symbol']) for s in spans),n)
        tif=ROOT/'figures'/f'Fig{n}.tif'
        with Image.open(tif) as im:
            check('figure_TIFF_contract',im.mode=='RGB' and min(im.info.get('dpi',(0,0)))>=599.5 and im.n_frames==1 and tif.stat().st_size<10000000,n)
    for path in ROOT.rglob('*.py'):
        if '__pycache__' in path.parts:continue
        ast.parse(path.read_text(encoding='utf-8-sig'),filename=str(path))
    check('all_python_AST',True,len(list(ROOT.rglob('*.py'))))
    cfr=json.loads((AUDIT/'cfr_translation_audit.json').read_text(encoding='utf-8'))
    check('CFR_audit_passed',not cfr['failed'] and cfr['passed']==cfr['checks'],cfr['passed'])
    for relative,digest in cfr['input_sha256'].items():
        check('CFR_saved_input_unchanged',sha(ROOT/relative)==digest,relative)
    for relative,digest in cfr['source_sha256'].items():
        check('CFR_model_source_unchanged',sha(ROOT/relative)==digest,relative)
    unchanged=0
    if args.baseline:
        baseline=args.baseline.resolve()
        check('citation_keys_and_occurrences_preserved',citations((baseline/'paper/main.tex').read_text())==used)
        for folder in ['research/inputs','research/outputs','research/coverage/inputs','research/coverage/outputs']:
            for original in (baseline/folder).rglob('*'):
                if not original.is_file() or '__pycache__' in original.parts:continue
                relative=original.relative_to(baseline)
                check('frozen_numerical_asset_unchanged',(ROOT/relative).is_file() and sha(original)==sha(ROOT/relative),relative.as_posix())
                unchanged+=1
        for relative in ['research/src/experiments.py','research/src/manuscript.tex']:
            if (baseline/relative).exists():check('retained_method_source',sha(ROOT/relative)==sha(baseline/relative),relative)
        old=(baseline/'paper/main.tex').read_text()
        def principal(s):
            section=s[s.index('\\label{tab:factorial}'):]
            return re.sub(r'\s+','',section.split('\\end{table}',1)[0])
        check('principal_result_table_unchanged',principal(old)==principal(text))
    failed=[x for x in checks if not x['passed']]
    result=dict(revision='v1.0.1',checks=len(checks),passed=len(checks)-len(failed),
                failed=failed,pages=page_count,citation_keys=len(used),
                original_numerical_files_compared=unchanged,
                complete_144_ODE_rerun=False,check_records=checks)
    (AUDIT/'revision_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='check_records'},indent=2))
    if failed:raise AssertionError(failed)


if __name__=='__main__':
    main()
