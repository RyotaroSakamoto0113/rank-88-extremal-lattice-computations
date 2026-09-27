from pathlib import Path
import argparse,subprocess,shutil,time,os,json,datetime
O=Path(__file__).resolve().parents[1]
a=argparse.ArgumentParser();a.add_argument('source');a.add_argument('name');args=a.parse_args()
R=O/'results'/args.name;R.mkdir(parents=True,exist_ok=False)
launcher=R/'launcher.gp';launcher.write_text('OUT='+json.dumps(str(R))+';\nROOT='+json.dumps(str(O))+';\n'+(O/'src'/args.source).read_text())
t=time.perf_counter()
with (R/'stdout').open('w') as so,(R/'stderr').open('w') as se:
 p=subprocess.Popen([shutil.which('gp'),'-fq',str(launcher)],cwd=O/'workspace',stdout=so,stderr=se)
 _,st,ru=os.wait4(p.pid,0);p.returncode=os.waitstatus_to_exitcode(st)
r={'returncode':p.returncode,'wall_seconds':time.perf_counter()-t,'user_cpu_seconds':ru.ru_utime,'system_cpu_seconds':ru.ru_stime,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
r['status']='verified' if not p.returncode and ('INITIAL78_VERIFIED' if args.source=='initial_check.gp' else 'COMPLETE') in (R/'stdout').read_text() and not any(line.lstrip().startswith('***') and 'Warning:' not in line for line in (R/'stderr').read_text().splitlines()) else 'failed'
(R/'timing.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r));print((R/'stdout').read_text());print((R/'stderr').read_text());raise SystemExit(r['status']!='verified')
