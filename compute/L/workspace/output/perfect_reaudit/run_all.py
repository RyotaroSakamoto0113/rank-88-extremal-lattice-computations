#!/usr/bin/env python3
"""Fresh exact reconstruction, deterministic seed search and perfectness proof."""
import argparse,datetime,hashlib,json,pathlib,resource,shutil,subprocess,sys,tempfile,time
H=pathlib.Path(__file__).resolve().parent;PROJECT=H.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(x,s):
 if not x:raise RuntimeError(s)
def main():
 p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=41);p.add_argument('--seed',type=int,default=20260913);p.add_argument('--gp',default='gp');p.add_argument('--cxx',default='clang++');p.add_argument('--output-root',type=pathlib.Path,help='write all new run files here and leave source pointers unchanged');a=p.parse_args()
 destination=a.output_root.resolve() if a.output_root else H/'runs';destination.mkdir(parents=True,exist_ok=True)
 run=pathlib.Path(tempfile.mkdtemp(prefix=datetime.datetime.now().strftime('%Y%m%d-%H%M%S-full-'),dir=destination));work=run/'workspace';out=work/'output/perfect_reaudit';(out/'input').mkdir(parents=True);(out/'src').mkdir();(out/'search').mkdir()
 names=['rank88_exact_data.gp','output/d3/benchmark/ambient.gp']
 for name in names:
  dest=work/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(PROJECT/name,dest)
 for rel in ['input/build_tensor_input.gp','input/verify_tensor_input.py','src/find_seeds.gp','src/split_rank.cpp','src/export_input.py','src/verify_split_certificate.py']:
  shutil.copy2(H/rel,out/rel)
 shutil.copy2(pathlib.Path(__file__),run/'run_all.py')
 info={'status':'running','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'trials':a.trials,'seed':a.seed,'run_directory':str(run),'stages':[],'initial_source_hashes':{str(p.relative_to(work)):sha(p)for p in work.rglob('*')if p.is_file()}}
 start=time.perf_counter();print('RUN_DIRECTORY',run,flush=True)
 def save():(run/'FULL_SUMMARY.json').write_text(json.dumps(info,indent=2)+'\n')
 def stage(name,cmd,marker=None):
  t=time.perf_counter();c=resource.getrusage(resource.RUSAGE_CHILDREN);row={'stage':name,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':list(map(str,cmd))};print('START',name,flush=True)
  with(run/(name+'.stdout')).open('w')as f,(run/(name+'.stderr')).open('w')as g:r=subprocess.run(list(map(str,cmd)),cwd=work,stdout=f,stderr=g)
  v=resource.getrusage(resource.RUSAGE_CHILDREN);row.update(returncode=r.returncode,wall_seconds=time.perf_counter()-t,user_cpu_seconds=v.ru_utime-c.ru_utime,system_cpu_seconds=v.ru_stime-c.ru_stime,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());info['stages'].append(row);save();need(r.returncode==0,name+' unsuccessful')
  errors=[x for x in(run/(name+'.stderr')).read_text().splitlines()if x.lstrip().startswith('***')and'Warning:'not in x];need(not errors,name+' GP error')
  if marker:need(marker in(run/(name+'.stdout')).read_text(),name+' missing completion marker')
  print('DONE',name,row['wall_seconds'],flush=True)
 try:
  stage('exact_input_reconstruction',[a.gp,'-fq',out/'input/build_tensor_input.gp'],'TENSOR_INPUT_VERIFIED')
  stage('portable_input_export',[sys.executable,out/'src/export_input.py',out/'input'])
  stage('independent_field_input_verification',[sys.executable,out/'input/verify_tensor_input.py'])
  need(json.loads((out/'input/independent_verification.json').read_text())['status']=='verified','field input verification')
  need(json.loads((out/'input/lattice.json').read_text())==json.loads((H/'input/lattice.json').read_text()),'different source lattice')
  stage('compile_split_rank',[a.cxx,'-O3','-std=c++17',out/'src/split_rank.cpp','-o',out/'src/split_rank'])
  launch='default(nbthreads,1);default(parisizemax,3000000000);default(parisize,256000000);\n'+f'OUT={json.dumps(str(out/"search"))};TRIALS={a.trials};RNGSEED={a.seed};\n'+f'read({json.dumps(str(out/"input/lattice.gp"))});read({json.dumps(str(out/"src/find_seeds.gp"))});\n'
  (out/'search/launch.gp').write_text(launch)
  stage('seed_search',[a.gp,'-fq',out/'search/launch.gp'],'SEED_SEARCH_COMPLETE')
  candidates=[json.loads(x)for x in(out/'search/candidates.jsonl').read_text().splitlines()]
  with(out/'search/candidates.txt').open('w')as f:
   f.write(str(len(candidates))+' 88\n');f.writelines(' '.join(map(str,v['x']))+'\n'for v in candidates)
  stage('split_rank',[out/'src/split_rank',out/'input/integer_input.txt',out/'search/candidates.txt',out/'certificate.json',461],'SPLIT_RANK_RESULT 3916')
  stage('independent_perfection_certificate',[sys.executable,out/'src/verify_split_certificate.py','--input',out/'input/lattice.json','--certificate',out/'certificate.json','--output',out/'verification'])
  cert=json.loads((out/'certificate.json').read_text());v=json.loads((out/'verification/verification.json').read_text());need(v['status']=='verified'and v['rational_span_dimension']==3916,'independent fullrank')
  info.update(status='verified',seed_count=len(cert['seeds']),candidate_count=len(candidates),rational_span_dimension=3916,modular_prime=461,certificate=str(out/'certificate.json'),certificate_sha256=sha(out/'certificate.json'),input_sha256=sha(out/'input/lattice.json'),independent_verification=str(out/'verification/verification.json'),external_minimum_8_prerequisite=True)
  if not a.output_root:(H/'LATEST_FULL.txt').write_text(str(run)+'\n')
 except BaseException as e:info.update(status='failed',error=str(e));raise
 finally:info.update(wall_seconds=time.perf_counter()-start,user_cpu_seconds=sum(s['user_cpu_seconds']for s in info['stages']),system_cpu_seconds=sum(s['system_cpu_seconds']for s in info['stages']),finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());save()
 print(json.dumps({k:v for k,v in info.items()if k not in ['stages','initial_source_hashes']},indent=2))
if __name__=='__main__':main()
