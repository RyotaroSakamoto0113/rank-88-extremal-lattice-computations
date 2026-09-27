"""Recompute the projective minima and save every final certificate."""
from pathlib import Path
import subprocess,sys,os,time,json,datetime,shutil
O=Path(__file__).resolve().parents[1]; W=O/'workspace'; A=W/'output/minima_reaudit'; D=W/'output/d4_projective_J'
(A/'runs').mkdir(exist_ok=True)
R=O/'results/core'; R.mkdir(parents=True,exist_ok=False)
summary={'status':'running','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stages':[]}
start=time.perf_counter()
def save(): (R/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
def run(name,command):
    command=list(map(str,command));t=time.perf_counter()
    print('START',name,flush=True)
    with (R/(name+'.stdout')).open('w') as so,(R/(name+'.stderr')).open('w') as se:
        p=subprocess.Popen(command,cwd=W,stdout=so,stderr=se)
        _,status,ru=os.wait4(p.pid,0);p.returncode=os.waitstatus_to_exitcode(status)
    row={'stage':name,'command':command,'returncode':p.returncode,'wall_seconds':time.perf_counter()-t,'user_cpu_seconds':ru.ru_utime,'system_cpu_seconds':ru.ru_stime}
    summary['stages'].append(row);save()
    if p.returncode:raise RuntimeError(name+' failed')
    print('DONE',name,round(row['wall_seconds'],3),flush=True)
def py(name,rel,*args):run(name,[sys.executable,A/rel,*args])
try:
    py('bootstrap','bootstrap_run.py')
    py('d1_d2_d3','d3_run.py')
    py('independent_projective_d3','theory/run_independent_d3.py')
    py('ordinary_sphere','ordinary_run.py')
    ordinary=Path((A/'LATEST_ORDINARY.txt').read_text().strip())/'output/rank4_audit'
    py('lookup','witnesses/ordinary_prune/build_lookup.py','--sphere',ordinary/'M_shell12.txt')
    reps=W/'output/rank4_d27_d48/runs/20260912-175210-2mzyf0pb/output/rank4_audit/M_shell12_orbit_reps.txt'
    py('upper_witnesses','witnesses/run_witnesses.py','--shell',ordinary/'M_shell12.txt','--orbits',reps,'--output',A/'witnesses/fresh_reproduction')
    py('initial_lines','initial4/run_initial4.py')
    py('initial_stabilizers','initial4/build_stab4.py')
    boot=Path((A/'LATEST_BOOTSTRAP.txt').read_text().strip())
    for src,dst in [(boot/'bootstrap/ambient.gp',D/'inputs/ambient.gp'),(boot/'bootstrap/generators.gp',D/'inputs/generators.gp'),
                    (boot/'exact_enum',D/'inputs/exact_enum'),(A/'initial4/FIRST4.gp',D/'inputs/FIRST4.gp'),
                    (A/'initial4/STAB4.gp',D/'inputs/STAB4.gp')]:
        shutil.copy2(src,dst);dst.chmod(0o755 if dst.name=='exact_enum' else 0o644)
    for name in ['ordinary_lookup.gp','shorter_ordinary.gp']:
        shutil.copy2(A/'witnesses/ordinary_prune'/name,D/'inputs/ordinary_prune'/name)
    if '--prerequisites-only' in sys.argv:
        summary.update(status='prerequisites_verified',wall_seconds=time.perf_counter()-start);save();print('CORE_PREREQUISITES_COMPLETE',flush=True);sys.exit(0)
    run('projective_d4_all78',[sys.executable,D/'run_full.py','--jobs','4'])
    d4=json.loads((D/'full_run/SUMMARY.json').read_text())
    assert d4['status']=='verified' and d4['exclusion_complete'] and d4['counts']['raw_rank4']==0
    summary.update(status='verified',d4_projective_lower=48,all78_verified=True,
        bootstrap=str(boot.relative_to(W)),d3=(A/'LATEST_D3.txt').read_text().strip(),
        independent_d3=(A/'theory/LATEST_INDEPENDENT_D3.txt').read_text().strip(),
        ordinary=(A/'LATEST_ORDINARY.txt').read_text().strip(),d4_counts=d4['counts'])
except SystemExit:
    raise
except BaseException as e:
    summary.update(status='failed',error=type(e).__name__+': '+str(e));raise
finally:
    summary.update(wall_seconds=time.perf_counter()-start,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());save()
print('CORE_COMPLETE',summary['status'],summary['wall_seconds'],flush=True)
