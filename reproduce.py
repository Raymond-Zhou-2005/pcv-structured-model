"""Reproduce the paper's numerical analyses and all main/extended figures.

Run from any working directory. Outputs stay within this repository.
The default pipeline uses retained closed-cohort replicates; --with-cohort
regenerates that optional 4,000-replicate audit with its recorded seed.
"""
from pathlib import Path
import argparse, json, shutil, subprocess, sys, time

ROOT=Path(__file__).resolve().parent
LOG=ROOT/'checks'/'logs'

def run(relative,*args):
    LOG.mkdir(parents=True,exist_ok=True)
    path=ROOT/relative
    label=relative.replace('/','_').removesuffix('.py')
    log=LOG/(label+'.log')
    command=[sys.executable,'-B',str(path),*map(str,args)]
    print('Running',relative,flush=True)
    with log.open('w',encoding='utf-8') as out:
        result=subprocess.run(command,cwd=ROOT,stdout=out,stderr=subprocess.STDOUT)
    if result.returncode:
        print(log.read_text(encoding='utf-8')[-6000:]);raise SystemExit(result.returncode)
    print('PASS',relative,flush=True)

def refresh_plot_data():
    dest=ROOT/'revision/data';dest.mkdir(exist_ok=True)
    for name in ['ctmc.csv','manufactured.csv','refinement.csv','aggregation.csv',
                 'aggregation_states.npz','aom.csv','mortality_brackets.csv']:
        shutil.copy2(ROOT/'research/outputs'/name,dest/name)
    for name in ['factorial.csv','sensitivity.csv']:
        shutil.copy2(ROOT/'research/inputs/Final/outputs'/name,dest/name)
    shutil.copy2(ROOT/'research/inputs/Final/outputs/runs/pcv13_D1W1G1/parameters.json',dest/'parameters.json')
    shutil.copy2(ROOT/'research/outputs/ctmc.csv',ROOT/'revision/inputs/ctmc.csv')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['all','analysis','plots','check'],default='all')
    parser.add_argument('--with-cohort',action='store_true')
    args=parser.parse_args();start=time.time()
    if args.stage in ['analysis','all']:
        run('research/inputs/Final/src/analyze.py')
        run('research/inputs/Final/src/validate.py')
        if args.with_cohort:
            run('research/inputs/Cohort/run_ode_microsim_reconciliation.py',
                '--output-dir',ROOT/'research/outputs/cohort','--replicates','1000','--seed','20261004')
        run('research/src/experiments.py')
        run('research/coverage/reproduce.py')
        run('research/sources/audit_sources.py')
        refresh_plot_data()
        run('revision/ctmc_variability.py')
    if args.stage in ['plots','all']:
        refresh_plot_data()
        run('revision/plot_figures.py')
        run('plots.py')
    run('checks/revision_20261103/cfr_translation_audit.py')
    run('verify.py')
    (ROOT/'checks/last_run.json').write_text(json.dumps(dict(
        stage=args.stage,with_cohort=args.with_cohort,elapsed_seconds=time.time()-start,
        python=sys.version,script=str(Path(__file__).relative_to(ROOT)),
        expected_scope='Deterministic scenario reproduction and conditional simulation checks; not empirical calibration or PSA'),indent=2),encoding='utf-8')

if __name__=='__main__':main()
