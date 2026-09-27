#!/usr/bin/env python3
"""Independently verify every completed D4 shard, with bounded concurrency.

At most two verifier processes run while production is active. Once the
production run finishes, the requested post-production concurrency is used.
Every verifier's child/parent CPU times and its concurrent wall interval are
preserved. The aggregate is successful only if all78 anchors verify.
"""
import argparse,concurrent.futures,datetime,hashlib,json,os,shutil,subprocess,sys,tempfile,time
from pathlib import Path

def need(ok,s):
    if not ok:raise ValueError(s)

def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def atomic(p,x):
    tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(p)

def read_json(p):
    try:return json.loads(p.read_text())
    except (FileNotFoundError,json.JSONDecodeError):return None

def main():
    src=Path(__file__).resolve().parent
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run',type=Path);p.add_argument('--jobs',type=int,default=4);p.add_argument('--while-production-jobs',type=int,default=2);p.add_argument('--gp',default='gp');a=p.parse_args()
    need(1<=a.while_production_jobs<=2 and a.jobs>=a.while_production_jobs,'valid concurrency limits');production=a.run.resolve();need((production/'anchors').is_dir(),'production anchors directory')
    audit=Path(tempfile.mkdtemp(prefix='independent-verification-',dir=production));(audit/'src').mkdir();(audit/'anchors').mkdir();gp=shutil.which(a.gp) or a.gp
    for name in ['verify_d4_run.py','verify_d4_shard.py','verify_d4_shard.gp','verify_d4_shard_README.md']:
        shutil.copy2(src/name,audit/'src'/name);(audit/'src'/name).chmod(0o444)
    lookup=production/'inputs/ordinary_prune/ordinary_lookup.gp';need(lookup.is_file(),'snapshotted ordinary lookup')
    started=time.perf_counter();summary={'status':'running','started_utc':utc(),'production_run':str(production),'verification_run':str(audit),'concurrency_while_production':a.while_production_jobs,'concurrency_after_production':a.jobs,'input_hashes':{str(x):digest(x) for x in [lookup,*sorted((audit/'src').iterdir())]},'anchors':{},'errors':[]}
    def save():
        records=list(summary['anchors'].values());good=[r for r in records if r.get('status')=='verified'];summary.update(elapsed_wall_including_waits_seconds=time.perf_counter()-started,verified_anchor_ids=sorted(r['anchor'] for r in good),sum_verification_wall_seconds=sum(r.get('total_wall_seconds',0) for r in good),sum_user_cpu_seconds=sum(r.get('gp_stage',{}).get('user_cpu_seconds',0)+r.get('python_user_cpu_seconds',0) for r in records),sum_system_cpu_seconds=sum(r.get('gp_stage',{}).get('system_cpu_seconds',0)+r.get('python_system_cpu_seconds',0) for r in records));summary['counts']=[sum(r['counts'][i] for r in good) for i in range(13)];atomic(audit/'SUMMARY.json',summary)
    def verify(anchor):
        shard=production/'anchors'/f'{anchor:03d}';out=audit/'anchors'/f'{anchor:03d}';stamp=utc();t=time.perf_counter();command=[sys.executable,str(audit/'src/verify_d4_shard.py'),str(shard),'--output',str(out),'--lookup',str(lookup),'--gp',gp]
        stdout_path=audit/'anchors'/f'{anchor:03d}.stdout';stderr_path=audit/'anchors'/f'{anchor:03d}.stderr'
        with stdout_path.open('w') as stdout,stderr_path.open('w') as stderr:
            child=subprocess.Popen(command,stdout=stdout,stderr=stderr)
            while True:
                try:pid,status,usage=os.wait4(child.pid,0);break
                except InterruptedError:pass
            child.returncode=os.waitstatus_to_exitcode(status)
        record=read_json(out/'verification.json') or {'status':'failed','error':'verification result absent'}
        if record.get('anchor')!=anchor:record.update(status='failed',error='verified anchor does not match dispatched anchor')
        record['anchor']=anchor;record['dispatch']={'started_utc':stamp,'finished_utc':utc(),'wall_seconds':time.perf_counter()-t,'returncode':child.returncode,'wait4_user_cpu_seconds':usage.ru_utime,'wait4_system_cpu_seconds':usage.ru_stime,'max_rss_native_units':usage.ru_maxrss,'stdout':str(stdout_path),'stderr':str(stderr_path),'command':command}
        if child.returncode!=0 or record.get('status')!='verified':record['status']='failed';record.setdefault('error','verifier did not complete successfully')
        atomic(out/'dispatcher_result.json',record) if out.is_dir() else atomic(audit/'anchors'/f'{anchor:03d}.failure.json',record)
        return record
    pending={};dispatched=set();terminal=False;failure=False;pool=concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs)
    print('VERIFICATION_RUN',audit,flush=True);(src.parent/'LATEST_D4_VERIFICATION.txt').write_text(str(audit)+'\n');save()
    try:
        while True:
            state=read_json(production/'D4_SUMMARY.json');terminal=state is not None and state.get('status')!='running';limit=a.jobs if terminal else a.while_production_jobs
            if not failure:
                for anchor in range(78):
                    if len(pending)>=limit:break
                    if anchor in dispatched:continue
                    shard_state=read_json(production/'anchors'/f'{anchor:03d}'/'ANCHOR_SUMMARY.json')
                    if shard_state and shard_state.get('status') in ('complete','counterexample'):
                        dispatched.add(anchor);pending[pool.submit(verify,anchor)]=anchor
            if not pending:
                if failure or len(dispatched)==78 or terminal:break
                time.sleep(1);continue
            done,_=concurrent.futures.wait(pending,timeout=2,return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                anchor=pending.pop(future)
                try:record=future.result()
                except BaseException as e:record={'anchor':anchor,'status':'failed','error':type(e).__name__+': '+str(e)}
                summary['anchors'][str(anchor)]=record
                if record['status']!='verified':failure=True;summary['errors'].append({'anchor':anchor,'error':record.get('error')})
                print('ANCHOR_VERIFICATION',anchor,record['status'],'wall',record.get('total_wall_seconds'),'counts',record.get('counts'),flush=True)
            save()
    except BaseException as e:
        failure=True;summary['errors'].append({'error':type(e).__name__+': '+str(e)})
    finally:
        pool.shutdown(wait=True)
        for future,anchor in pending.items():
            try:record=future.result()
            except BaseException as e:record={'anchor':anchor,'status':'failed','error':str(e)}
            summary['anchors'][str(anchor)]=record
            if record['status']!='verified':failure=True;summary['errors'].append({'anchor':anchor,'error':record.get('error')})
    complete=len(summary['anchors'])==78 and all(r['status']=='verified' for r in summary['anchors'].values());summary.update(status='failed' if failure else 'verified' if complete else 'incomplete',finished_utc=utc(),all78_verified=complete,not_dispatched_anchor_ids=sorted(set(range(78))-dispatched));save();print('VERIFICATION_FINISHED',summary['status'],'elapsed',summary['elapsed_wall_including_waits_seconds'],'user_cpu',summary['sum_user_cpu_seconds'],'system_cpu',summary['sum_system_cpu_seconds'],flush=True)
    if not complete:raise SystemExit(1)

if __name__=='__main__':main()
