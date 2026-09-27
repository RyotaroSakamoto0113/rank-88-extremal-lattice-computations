#!/usr/bin/env python3
import argparse,datetime,hashlib,json,pathlib,resource,shutil,subprocess,tempfile,time
H=pathlib.Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=16);p.add_argument('--seed',type=int,default=20260913);a=p.parse_args()
 out=pathlib.Path(tempfile.mkdtemp(prefix=datetime.datetime.now().strftime('%Y%m%d-%H%M%S-seeds-'),dir=H/'runs'))
 for src in [H/'src/find_seeds.gp',H/'input/lattice.gp',H/'input/lattice.json',H/'input/integer_input.txt']:
  shutil.copy2(src,out/src.name)
 launch='default(nbthreads,1);default(parisizemax,3000000000);default(parisize,256000000);\n'+f'OUT={json.dumps(str(out))};TRIALS={a.trials};RNGSEED={a.seed};\n'+f'read({json.dumps(str(out/"lattice.gp"))});read({json.dumps(str(out/"find_seeds.gp"))});\n'
 (out/'launch.gp').write_text(launch);info={'status':'running','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'trials':a.trials,'seed':a.seed,'input_hashes':{x.name:sha(x) for x in out.iterdir() if x.is_file()}}
 print('RUN_DIRECTORY',out,flush=True);start=time.perf_counter();c=resource.getrusage(resource.RUSAGE_CHILDREN)
 with(out/'gp.stdout').open('w')as f,(out/'gp.stderr').open('w')as g:r=subprocess.run(['gp','-fq',str(out/'launch.gp')],stdout=f,stderr=g)
 v=resource.getrusage(resource.RUSAGE_CHILDREN);info.update(returncode=r.returncode,wall_seconds=time.perf_counter()-start,user_cpu_seconds=v.ru_utime-c.ru_utime,system_cpu_seconds=v.ru_stime-c.ru_stime,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 errors=[x for x in(out/'gp.stderr').read_text().splitlines() if x.lstrip().startswith('***') and 'Warning:'not in x]
 if r.returncode or errors or 'SEED_SEARCH_COMPLETE' not in(out/'gp.stdout').read_text():
  info.update(status='failed',errors=errors);(out/'SEARCH_SUMMARY.json').write_text(json.dumps(info,indent=2)+'\n');raise RuntimeError('seed search failed; see '+str(out))
 rows=[json.loads(x) for x in(out/'candidates.jsonl').read_text().splitlines()];info.update(status='complete',candidate_count=len(rows))
 with(out/'candidates.txt').open('w')as f:
  f.write(f'{len(rows)} 88\n');f.writelines(' '.join(map(str,x['x']))+'\n'for x in rows)
 (out/'SEARCH_SUMMARY.json').write_text(json.dumps(info,indent=2)+'\n');(H/'LATEST_SEARCH.txt').write_text(str(out)+'\n');print(json.dumps(info,indent=2))
if __name__=='__main__':main()
