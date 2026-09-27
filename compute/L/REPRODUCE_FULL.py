#!/usr/bin/env python3
"""Reproduce in a fresh directory, using only files inside this supplement."""
from pathlib import Path
import argparse,sys,subprocess,shutil,json,time,os,shlex
SOURCE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);ap.add_argument('--smoke',action='store_true',help='Portable prerequisites, h=3, perfection and group checks; omit long full enumerations and field class groups.');args=ap.parse_args()
O=args.output.resolve();O.mkdir(parents=True,exist_ok=False)
gp=shutil.which('gp');cxx=shutil.which(os.environ.get('CXX','clang++')) or shutil.which('g++');assert gp and cxx,'PARI/GP and a C++17 compiler are required'
for d in ['src','references']:shutil.copytree(SOURCE/d,O/d,ignore=shutil.ignore_patterns('__pycache__'))
files=json.loads((SOURCE/'SOURCE_MANIFEST.json').read_text())
for name in files:
 p=O/'workspace'/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(SOURCE/'workspace'/name,p)
 # Only executable discovery is changed; mathematical source is preserved.
 if p.suffix=='.py':p.write_text(p.read_text().replace('/opt/homebrew/bin/gp',gp))
(O/'results').mkdir();(O/'bin').mkdir();shim=O/'bin/clang++';shim.write_text('#!/bin/sh\nexec '+shlex.quote(cxx)+' "$@"\n');shim.chmod(0o755)
env=dict(os.environ,PATH=str(O/'bin')+os.pathsep+os.environ.get('PATH',''))
timings=[];start=time.perf_counter()
def run(name,cmd):
 print('START',name,flush=True);t=time.perf_counter();logs=O/'results/driver';logs.mkdir(exist_ok=True)
 with (logs/(name+'.stdout')).open('w') as so,(logs/(name+'.stderr')).open('w') as se:
  p=subprocess.Popen(list(map(str,cmd)),cwd=O,env=env,stdout=so,stderr=se);_,st,ru=os.wait4(p.pid,0);p.returncode=os.waitstatus_to_exitcode(st)
 r={'stage':name,'returncode':p.returncode,'wall_seconds':time.perf_counter()-t,'user_cpu_seconds':ru.ru_utime,'system_cpu_seconds':ru.ru_stime};timings.append(r);(O/'REPRODUCTION_TIMINGS.json').write_text(json.dumps(timings,indent=2)+'\n');print('DONE',name,r['returncode'],round(r['wall_seconds'],3),flush=True)
 if p.returncode:raise RuntimeError(name+' failed; see results/driver/'+name+'.stderr')
def py(name,p,*a):run(name,[sys.executable,O/p,*a])
def gpstage(source,name):py(name,'src/run_gp.py',source,name)
def compile(name,source,gmp=False):
 cmd=[cxx,'-O3','-std=c++17']
 if gmp and Path('/opt/homebrew/include').exists():cmd+=['-I/opt/homebrew/include','-L/opt/homebrew/lib']
 cmd+=[O/'src'/source]
 if gmp:cmd+=['-lgmpxx','-lgmp']
 cmd+=['-o',O/'results'/name];run('compile_'+name,cmd)
try:
 py('core','src/run_core.py',*(['--prerequisites-only'] if args.smoke else []))
 py('h3','workspace/output/rank4_d27_d48/run_h3.py','--run-dir',O/'results/h3')
 py('perfect','workspace/output/perfect_reaudit/reverify.py','--output',O/'results/perfect')
 py('additional','src/verify_additional.py')
 A=O/'workspace/output/minima_reaudit';b=Path((A/'LATEST_BOOTSTRAP.txt').read_text().strip());uv=json.loads((A/'witnesses/fresh_reproduction/witnesses.jsonl').read_text().splitlines()[1])['coordinates'];(O/'results/additional/binary_input.gp').write_text('UV='+str(uv)+';\n')
 gpstage('binary_witness.gp','binary_witness');gpstage('initial_check.gp','initial_check')
 compile('group_structure','group_structure.cpp');(O/'results/group').mkdir()
 run('group_structure',[O/'results/group_structure',A/'witnesses/orbits/group_matrices.txt',b/'bootstrap/cusp1_integer_data.txt',O/'workspace/output/perfect_reaudit/input/J_M23.txt',O/'results/group'])
 py('group_action','src/verify_group_action.py')
 gpstage('positivity.gp','positivity');gpstage('barnes.gp','barnes');gpstage('small_ideals.gp','small_ideals')
 if not args.smoke:
  gpstage('field_arithmetic.gp','fields');gpstage('nebe_class_number.gp','nebe_class_number')
  gpstage('orbit_burnside.gp','orbit_table');compile('fast_histogram','fast_histogram.cpp',True);compile('exact_histogram','exact_histogram.cpp',True)
  lines=(O/'results/orbit_table/exact_forms.txt').read_text().splitlines();forms=[];i=0
  while i<len(lines):n=int(lines[i].split()[0]);forms.append('\n'.join(lines[i:i+n+1])+'\n');i+=n+1
  small=''.join(f for f in forms if int(f.split()[0])<22);(O/'results/orbit_table/nonidentity_forms.txt').write_text(small)
  run('fast_all',[O/'results/fast_histogram','--batch',O/'results/orbit_table/exact_forms.txt']);(O/'results/fast_orbit_all').mkdir();shutil.copy2(O/'results/driver/fast_all.stdout',O/'results/fast_orbit_all/stdout')
  run('fast_small',[O/'results/fast_histogram','--batch',O/'results/orbit_table/nonidentity_forms.txt']);shutil.copy2(O/'results/driver/fast_small.stdout',O/'results/fast_nonidentity.jsonl')
  run('rational_small',[O/'results/exact_histogram','--batch',O/'results/orbit_table/nonidentity_forms.txt','60']);(O/'results/orbit_exact_remaining').mkdir();shutil.copy2(O/'results/driver/rational_small.stdout',O/'results/orbit_exact_remaining/nonidentity.stdout')
  py('fast_regression','src/test_fast_histogram.py');py('orbit_verification','src/verify_orbit_table.py')
  gpstage('verify_catalogue.gp','catalogue');run('catalogue_exact',[O/'results/fast_histogram','--batch',O/'results/catalogue/exact_forms.txt']);shutil.copy2(O/'results/driver/catalogue_exact.stdout',O/'results/catalogue/exact_histograms.jsonl')
  py('bind_all','src/bind_certificates.py')
 report={'status':'smoke_verified' if args.smoke else 'fully_reproduced','output':str(O),'wall_seconds':time.perf_counter()-start,'stages':timings}
except BaseException as exc:
 report={'status':'failed','error':repr(exc),'wall_seconds':time.perf_counter()-start,'stages':timings};(O/'REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n');raise
(O/'REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n');print(report['status'],O,flush=True)
