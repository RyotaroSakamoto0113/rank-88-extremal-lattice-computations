#!/usr/bin/env python3
"""Rebuild the exact tensor input and independently verify it, saving timings."""
import datetime,hashlib,json,resource,subprocess,sys,time
from pathlib import Path

def read_matrix(p):
    t=list(map(int,p.read_text().split()));n,m=t[:2]
    if len(t)!=2+n*m:raise ValueError('matrix size mismatch')
    return[t[2+i*m:2+(i+1)*m] for i in range(n)]

def main():
    out=Path(__file__).resolve().parent;project=out.parents[2];summary={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stages':[]};start=time.perf_counter()
    def run(name,cmd):
        t=time.perf_counter();c=resource.getrusage(resource.RUSAGE_CHILDREN);stamp=datetime.datetime.now(datetime.timezone.utc).isoformat();r=subprocess.run(cmd,cwd=project,capture_output=True,text=True);d=resource.getrusage(resource.RUSAGE_CHILDREN);(out/(name+'.stdout')).write_text(r.stdout);(out/(name+'.stderr')).write_text(r.stderr);summary['stages'].append({'name':name,'started_utc':stamp,'wall_seconds':time.perf_counter()-t,'user_cpu_seconds':d.ru_utime-c.ru_utime,'system_cpu_seconds':d.ru_stime-c.ru_stime,'returncode':r.returncode});(out/'run_summary.json').write_text(json.dumps(summary,indent=2)+'\n');r.check_returncode();return r
    r=run('rebuild',['gp','-q','-f',str(out/'build_tensor_input.gp')])
    if 'TENSOR_INPUT_VERIFIED' not in r.stdout:raise ValueError('GP verification did not complete')
    x={'G':read_matrix(out/'G88.txt'),'M5':read_matrix(out/'M5.txt'),'M23':read_matrix(out/'M23.txt')};(out/'lattice.json').write_text(json.dumps(x,separators=(',',':'))+'\n');(out/'lattice.gp').write_text(''.join(name+'=['+';'.join(','.join(map(str,row)) for row in A)+'];\n' for name,A in x.items()));(out/'integer_input.txt').write_text('88\n'+'\n'.join(' '.join(map(str,row)) for A in x.values() for row in A)+'\n')
    run('independent_python',[sys.executable,str(out/'verify_tensor_input.py')]);check=json.loads((out/'independent_verification.json').read_text())
    if check['status']!='verified':raise ValueError('independent verification failed')
    summary.update(status='verified',total_wall_seconds=time.perf_counter()-start,lattice_sha256=hashlib.sha256((out/'lattice.json').read_bytes()).hexdigest());(out/'run_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
