"""All 78 anchors, all final classes, exact counts and independent algebra audit."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import datetime,importlib.util,json,sys,time,hashlib,shutil,traceback,argparse
O=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('d4_engine',O/'src/d4_run.py')
engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
def atomic(p,x):
    q=p.with_name(p.name+'.tmp');q.write_text(json.dumps(x,indent=2)+'\n');q.replace(p)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--jobs',type=int,default=4);args=ap.parse_args()
    run=O/'full_run';run.mkdir(exist_ok=False)
    shutil.copytree(O/'inputs',run/'inputs');shutil.copytree(O/'src',run/'src')
    files={str(f.relative_to(run)):engine.digest(f) for d in ['inputs','src'] for f in (run/d).rglob('*') if f.is_file()}
    atomic(run/'INPUTS_SHA256.json',files)
    runner=engine.Runner(run,shutil.which('gp'))
    t=time.perf_counter();summary={'status':'running','started_utc':engine.utcnow(),'bounds_by_final_class':[47,47,47],
       'ternary_hermite_cube_bound':46,'jobs':args.jobs,'anchors':{},'input_hashes':files}
    def save():
        summary['wall_seconds']=time.perf_counter()-t
        good=[r for r in summary['anchors'].values() if r['status']=='verified']
        summary['verified_anchors']=sorted(r['anchor_id'] for r in good)
        summary['counts']={name:sum(r['counts'][name] for r in good) for name in engine.COUNT_NAMES}
        summary['user_cpu_seconds']=sum(r.get('total_user_cpu_seconds',0) for r in summary['anchors'].values())
        summary['system_cpu_seconds']=sum(r.get('total_system_cpu_seconds',0) for r in summary['anchors'].values())
        atomic(run/'SUMMARY.json',summary)
    def work(anchor):
        record=runner.anchor(anchor)
        print('SEARCH',anchor,record['status'],'wall',round(record['wall_seconds'],3),'counts',record.get('counts'),flush=True)
        if record['status'] not in ('complete','counterexample'):
            return record
        shard=run/'anchors'/f'{anchor:03d}'
        command=[sys.executable,run/'src/verify_d4_shard.py',shard,'--output',shard/'independent_certificate',
                 '--lookup',run/'inputs/ordinary_prune/ordinary_lookup.gp','--gp',shutil.which('gp')]
        stage=runner.process(shard,'independent_algebra_audit',command)
        record['stages'].append(stage)
        try:
            proof=json.loads((shard/'independent_certificate/verification.json').read_text())
            assert stage['returncode']==0 and proof['status']=='verified'
            assert proof['anchor']==anchor and proof['counts']==record['raw_counts_vector'][:13]
            record['status']='verified';record['independent_algebra_verification']=proof
        except BaseException as exc:
            record.update(status='failed',error=str(exc))
        record['total_user_cpu_seconds']=sum(s['user_cpu_seconds'] for s in record['stages'])
        record['total_system_cpu_seconds']=sum(s['system_cpu_seconds'] for s in record['stages'])
        record['total_wall_seconds']=sum(s['wall_seconds'] for s in record['stages'])
        atomic(shard/'FINAL_ANCHOR.json',record)
        return record
    print('RUN',run,flush=True);save()
    # Bound concurrency including each shard's independent audit.
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures={pool.submit(work,n):n for n in range(78)}
        for f in as_completed(futures):
            n=futures[f]
            try:r=f.result()
            except BaseException as e:r={'anchor_id':n,'status':'failed','error':traceback.format_exc()}
            summary['anchors'][str(n)]=r;save()
            print('VERIFIED',n,r['status'],'done',len(summary['verified_anchors']),'/78','elapsed',round(summary['wall_seconds'],3),flush=True)
    complete=len(summary['verified_anchors'])==78
    summary.update(status='verified' if complete else 'failed',finished_utc=engine.utcnow(),
                   all78_verified=complete,exclusion_complete=complete and summary['counts']['raw_rank4']==0)
    save();print('FULL_RUN_FINISHED',summary['status'],'exclusion_complete',summary['exclusion_complete'],'wall',summary['wall_seconds'],flush=True)
    if not complete:raise SystemExit(1)
if __name__=='__main__':main()
