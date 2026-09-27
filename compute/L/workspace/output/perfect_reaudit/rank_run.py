#!/usr/bin/env python3
import argparse,datetime,hashlib,json,pathlib,resource,shutil,subprocess,tempfile,time
H=pathlib.Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--search',type=pathlib.Path);p.add_argument('--prime',type=int,default=461);a=p.parse_args();search=a.search or pathlib.Path((H/'LATEST_SEARCH.txt').read_text().strip())
 out=pathlib.Path(tempfile.mkdtemp(prefix=datetime.datetime.now().strftime('%Y%m%d-%H%M%S-split-'),dir=H/'runs'))
 for src in [H/'src/split_rank.cpp',search/'integer_input.txt',search/'candidates.txt',search/'lattice.json']:shutil.copy2(src,out/src.name)
 summary={'status':'running','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'search_directory':str(search),'prime':a.prime,'input_hashes':{f.name:sha(f)for f in out.iterdir()},'stages':[]};start=time.perf_counter()
 def run(name,cmd):
  t=time.perf_counter();c=resource.getrusage(resource.RUSAGE_CHILDREN);s={'stage':name,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':list(map(str,cmd))}
  with(out/(name+'.stdout')).open('w')as f,(out/(name+'.stderr')).open('w')as g:r=subprocess.run(list(map(str,cmd)),stdout=f,stderr=g)
  v=resource.getrusage(resource.RUSAGE_CHILDREN);s.update(returncode=r.returncode,wall_seconds=time.perf_counter()-t,user_cpu_seconds=v.ru_utime-c.ru_utime,system_cpu_seconds=v.ru_stime-c.ru_stime,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());summary['stages'].append(s);return r.returncode
 print('RUN_DIRECTORY',out,flush=True)
 try:
  if run('compile',['clang++','-O3','-std=c++17',out/'split_rank.cpp','-o',out/'split_rank']):raise RuntimeError('compilefailed')
  rc=run('split_rank',[out/'split_rank',out/'integer_input.txt',out/'candidates.txt',out/'certificate.json',a.prime]);cert=json.loads((out/'certificate.json').read_text());summary.update(status='complete'if rc==0 and cert['rank']==3916 else'partial',rank=cert['rank'],seed_count=len(cert['seeds']),processed_candidates=cert['timing']['processed_candidates'],certificate_sha256=sha(out/'certificate.json'))
  if summary['status']=='complete':(H/'LATEST_CERTIFICATE.txt').write_text(str(out)+'\n')
  if rc not in [0,2]:raise RuntimeError('rankfailed')
 except BaseException as e:summary.update(status='failed',error=str(e));raise
 finally:summary.update(wall_seconds=time.perf_counter()-start,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());(out/'RANK_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
