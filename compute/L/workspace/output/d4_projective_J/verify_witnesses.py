from pathlib import Path
import json,os,subprocess,datetime,time
from fractions import Fraction
O=Path(__file__).resolve().parent
WOUT=O/'witnesses/all_classes'
WOUT.mkdir(exist_ok=True,parents=True)
def transpose(A):return list(map(list,zip(*A)))
def mul(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def det(A):
    a=[list(map(Fraction,row)) for row in A];v=Fraction(1)
    for i in range(len(a)):
        k=next((k for k in range(i,len(a)) if a[k][i]),None)
        if k is None:return 0
        if k!=i:a[i],a[k]=a[k],a[i];v=-v
        q=a[i][i];v*=q
        for k in range(i+1,len(a)):
            b=a[k][i]/q
            for j in range(i+1,len(a)):a[k][j]-=b*a[i][j]
            a[k][i]=0
    return v
records=[];start=time.perf_counter()
for anchor in [0,1,'free']:
    probe=O/'probe/free' if anchor=='free' else next(x for x in sorted((O/'probe').glob(f'{anchor:03d}-*')) if (x/'witnesses.gp').stat().st_size)
    p=WOUT/f'anchor{anchor}'
    launcher=p.with_suffix('.gp')
    launcher.write_text('\n'.join(['default(parisizemax,3000000000);',
      'INPUTS='+json.dumps(str(O/'inputs'))+';',
      'WITNESS_FILE='+json.dumps(str(probe/'witnesses.gp'))+';',
      'OUT_JSON='+json.dumps(str(p.with_suffix('.json')))+';',
      'OUT_GP='+json.dumps(str(p.with_suffix('.pseudobasis.gp')))+';',
      'read('+json.dumps(str(O/'src/export_witness.gp'))+');'])+'\n')
    t=time.perf_counter()
    with p.with_suffix('.stdout').open('w') as so,p.with_suffix('.stderr').open('w') as se:
        proc=subprocess.Popen(['/opt/homebrew/bin/gp','-fq',str(launcher)],stdout=so,stderr=se)
        pid,st,ru=os.wait4(proc.pid,0);proc.returncode=os.waitstatus_to_exitcode(st)
    record={'anchor':anchor,'wall_seconds':time.perf_counter()-t,'user_cpu_seconds':ru.ru_utime,'system_cpu_seconds':ru.ru_stime,'returncode':proc.returncode}
    records.append(record)
    assert proc.returncode==0 and 'WITNESS_VERIFIED' in p.with_suffix('.stdout').read_text(),p.with_suffix('.stderr').read_text()
    r=json.loads(p.with_suffix('.json').read_text());B=r['basis'];G=r['ambient_trace_gram'];W=r['ambient_omega_action'];A=r['omega_action'];Q=r['trace_gram']
    assert mul(B,A)==mul(W,B) and mul(mul(transpose(B),G),B)==Q
    assert det(Q)==23**4*48**2
    assert r['steinitz_class_mod3']=={0:2,1:1,'free':0}[anchor]
    record.update(status='verified',trace_determinant=int(det(Q)),hermitian_determinant=48,steinitz_class_mod3=r['steinitz_class_mod3'],python_integer_and_rational_recheck=True)
    print(record,flush=True)
summary={'status':'verified','witnesses':records,'wall_seconds':time.perf_counter()-start,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(WOUT/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
