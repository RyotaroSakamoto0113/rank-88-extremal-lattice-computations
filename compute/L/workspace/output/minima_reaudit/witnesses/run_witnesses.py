#!/usr/bin/env python3
"""Compile, rediscover and independently verify the four upper witnesses."""
import argparse,datetime,hashlib,json,platform,resource,shutil,subprocess,sys,time
from pathlib import Path

def main():
    root=Path(__file__).resolve().parent
    legacy=root.parents[2]/'output/rank4_d27_d48/runs/20260912-175210-2mzyf0pb/output/rank4_audit'
    p=argparse.ArgumentParser();p.add_argument('--ambient',type=Path,default=root/'M_integer_data.txt');p.add_argument('--shell',type=Path,default=legacy/'M_shell12.txt');p.add_argument('--orbits',type=Path,default=legacy/'M_shell12_orbit_reps.txt');p.add_argument('--output',type=Path,default=root/'reproduction');a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    summary={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'platform':platform.platform(),'python':sys.version,'stages':[],'inputs':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [a.ambient,a.shell,a.orbits,root/'find_witnesses.cpp',root/'verify_witnesses.py']}}
    t0=time.perf_counter();c0=resource.getrusage(resource.RUSAGE_CHILDREN)
    def run(stage,args):
        started=datetime.datetime.now(datetime.timezone.utc).isoformat();t=time.perf_counter();c=resource.getrusage(resource.RUSAGE_CHILDREN)
        r=subprocess.run(args,capture_output=True,text=True);end=time.perf_counter();d=resource.getrusage(resource.RUSAGE_CHILDREN)
        (a.output/(stage+'.stdout')).write_text(r.stdout);(a.output/(stage+'.stderr')).write_text(r.stderr)
        summary['stages'].append({'stage':stage,'started_utc':started,'wall_seconds':end-t,'user_cpu_seconds':d.ru_utime-c.ru_utime,'system_cpu_seconds':d.ru_stime-c.ru_stime,'exit_code':r.returncode,'command':list(map(str,args))})
        (a.output/'timing.json').write_text(json.dumps(summary,indent=2)+'\n');r.check_returncode();return r
    exe=a.output/'find_witnesses';run('compile',[shutil.which('clang++') or shutil.which('g++'),'-O3','-std=c++17',str(root/'find_witnesses.cpp'),'-o',str(exe)])
    r=run('search',[str(exe),str(a.ambient),str(a.shell),str(a.orbits)])
    witnesses=a.output/'witnesses.jsonl';witnesses.write_text(r.stdout)
    run('independent_verification',[sys.executable,str(root/'verify_witnesses.py'),'--input',str(a.ambient),'--witnesses',str(witnesses),'--output',str(a.output/'verification.json')])
    d=resource.getrusage(resource.RUSAGE_CHILDREN);summary.update(status='verified',total_wall_seconds=time.perf_counter()-t0,total_child_user_cpu_seconds=d.ru_utime-c0.ru_utime,total_child_system_cpu_seconds=d.ru_stime-c0.ru_stime)
    (a.output/'timing.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
