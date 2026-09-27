#!/usr/bin/env python3
"""Regenerate the improved h=3 certificate; time every stage in a fresh directory."""
import argparse,datetime,hashlib,json,os,platform,re,resource,shutil,subprocess,sys,tempfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[1]
LEGACY=PROJECT/'output/rank4_audit'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(p):return [json.loads(x)for x in Path(p).read_text().splitlines()if x.strip()]
def flat(rr,p):Path(p).write_text(''.join(str(t['id'])+' '+' '.join(str(v)for r in t['H']for z in r for v in z)+'\n'for t in rr))
def gpm(H):return '['+';'.join(','.join(f'({a})+({b})*w'for a,b in row)for row in H)+']'
def need(b,s):
 if not b:raise RuntimeError(s)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run-dir',type=Path);ap.add_argument('--cached-ambient',action='store_true');ap.add_argument('--prepare-only',action='store_true');args=ap.parse_args()
 if args.run_dir:root=args.run_dir.resolve();root.mkdir(parents=True,exist_ok=False)
 else:
  (HERE/'runs').mkdir(exist_ok=True);root=Path(tempfile.mkdtemp(prefix=datetime.datetime.now().strftime('%Y%m%d-%H%M%S-'),dir=HERE/'runs'))
 out=root/'output/rank4_audit';out.mkdir(parents=True);logs=out/'logs';logs.mkdir()
 for f in LEGACY.iterdir():
  if f.is_file() and f.suffix in ('.cpp','.py','.gp'):shutil.copy2(f,out/f.name)
 for f in (HERE/'src').iterdir():
  if f.is_file() and f.suffix in ('.cpp','.py','.gp','.json'):shutil.copy2(f,out/f.name)
 shutil.copy2(HERE/'abstract_h3.cpp',out/'abstract_h3.cpp')
 deps=['rank88_exact_data.gp','output/d3/benchmark/ambient.gp','output/d3/benchmark/exact_enum.cpp',
       'output/direct_extremality/countermodel_data.gp','output/direct_extremality/basis_search_target_result.gp','output/d4/feasibility/free_counterexample.gp']
 for rel in deps:
  f=root/rel;f.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(PROJECT/rel,f)
 gp=shutil.which('gp') or '/opt/homebrew/bin/gp';cc=shutil.which('clang++') or 'g++'
 env=os.environ.copy();env['AUDIT_GP']=gp;env['PYTHONUNBUFFERED']='1'
 summary={'status':'prepared','run_directory':str(root),'platform':platform.platform(),'python':sys.version,'gp':gp,'compiler':cc,
  'bounds':{'d1_free':6,'d2_free':16,'d3_free':27,'d4_free':48,'provenance':'Provided as established by the user; the new minima are not reproved in this run.'},
  'h':3,'cached_ambient':args.cached_ambient,'driver_sha256':sha(__file__),'stages':[]}
 (out/'INPUT_MANIFEST.json').write_text(json.dumps({str(f.relative_to(root)):sha(f)for f in root.rglob('*')if f.is_file()},indent=2)+'\n')
 def save(): (out/'FULL_RUN_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
 save();print('RUN_DIRECTORY',root,flush=True)
 if args.prepare_only:return
 started=time.perf_counter();usage0=resource.getrusage(resource.RUSAGE_CHILDREN)
 def run(label,cmd,destination=None):
  print('START',label,flush=True);ts=time.perf_counter();before=resource.getrusage(resource.RUSAGE_CHILDREN)
  dest=Path(destination)if destination else logs/(label+'.stdout');err=logs/(label+'.stderr')
  with dest.open('w')as a,err.open('w')as b:r=subprocess.run(list(map(str,cmd)),cwd=root,env=env,stdout=a,stderr=b)
  after=resource.getrusage(resource.RUSAGE_CHILDREN)
  timing={'stage':label,'wall_seconds':time.perf_counter()-ts,'user_cpu_seconds':after.ru_utime-before.ru_utime,'system_cpu_seconds':after.ru_stime-before.ru_stime,
    'returncode':r.returncode,'stdout':str(dest.relative_to(root)),'stderr':str(err.relative_to(root))}
  summary['stages'].append(timing);save();need(r.returncode==0,f'{label} failed: {err}')
  if str(cmd[0])==gp:
   need(not any(line.lstrip().startswith('***')and'Warning:'not in line for line in err.read_text().splitlines()),f'PARI error in {label}')
  print('DONE',label,f"{timing['wall_seconds']:.3f}s",flush=True)
 def py(name,*a):run(Path(name).stem,[sys.executable,out/name,*a])
 def gpr(name):run(Path(name).stem,[gp,'-fq',out/name])
 flags=['-O3','-std=c++17'];gmp=[]
 for prefix in ('/opt/homebrew','/usr/local'):
  if Path(prefix+'/include/gmpxx.h').exists():gmp=['-I'+prefix+'/include','-L'+prefix+'/lib'];break
 gmp+=['-lgmpxx','-lgmp']
 try:
  for name,uses in [('abstract_h3',False),('check_M_shell_orbits',False),('check_binary_absence',False),('independent_embedding_sieve',True),('independent_simultaneous',True),('exact_affine_enum',True)]:
   run('compile_'+name,[cc,*flags,out/(name+'.cpp'),'-o',out/name,*(gmp if uses else [])])
  exact=root/'output/d3/benchmark/exact_enum';run('compile_exact_enum',[cc,*flags,root/'output/d3/benchmark/exact_enum.cpp','-o',exact,*gmp])
  gpr('verify_ambient_input.gp');need('AMBIENT_INPUT_VERIFIED'in(logs/'verify_ambient_input.stdout').read_text(),'ambient identity marker missing')
  (out/'generator_check.gp').write_text('read("output/d3/benchmark/ambient.gp");\nread("output/rank4_audit/ambient_qfauto.gp");\nfor(i=1,#AUT[2],A=AUT[2][i];if(denominator(A)!=1 || A~*G*A!=G || A*W!=W*A || abs(matdet(A))!=1,error("invalid generator")));\nprint("GENERATORS_VERIFIED");\nquit;\n')
  gpr('generator_check.gp');gpr('independent_chart_audit.gp');gpr('export_charts_json.gp')
  if args.cached_ambient:
   cached=['M_integer_data.txt','M_shell12.txt','M_shell12_count_input.txt','M_shell12_exact_count.json']
   summary['ambient_cache_hashes']={f:sha(LEGACY/f)for f in cached}
   for name in cached:shutil.copy2(LEGACY/name,out/name)
   need(json.loads((out/'M_shell12_exact_count.json').read_text())['complete'],'incomplete cached count')
  else:
   gpr('independent_M_shell.gp')
   for name in ['M_integer_data.txt','M_shell12.txt','M_shell12_count_input.txt']:
    f=out/name;f.write_text(re.sub(r'[\[\],;]',' ',f.read_text()))
   run('M_exact_sphere_count',[exact,out/'M_shell12_count_input.txt','60'],out/'M_shell12_exact_count.json')
  ec=json.loads((out/'M_shell12_exact_count.json').read_text());need(ec['complete']and ec['signed_vectors']==567226,'M sphere count')
  run('M_shell_orbits',[out/'check_M_shell_orbits'],out/'M_shell12_orbit_audit.json')
  orb=json.loads((out/'M_shell12_orbit_audit.json').read_text());need(orb['unique_valid']and orb['generator_closed']and sum(r['orbits']for r in orb['shells'])==78,'M sphere/orbit audit')
  need(2*orb['pair_count']==orb['signed_count']==ec['signed_vectors'],'M enumerated sphere/count equality')
  need(all(6<=r['q']<=12 and r['q']!=7 for r in orb['shells']),'M minimum and norm 7 absence')
  # The affine helper reads precisely the two integer matrices, followed by nothing else.
  v=(out/'M_integer_data.txt').read_text().split();need(int(v[0])==22,'M header');(out/'h2_M_input.txt').write_text('22 18\n'+' '.join(v[1:1+2*22*22])+'\n')
  run('abstract_h3',[out/'abstract_h3',out/'independent_h3','27','48','1'])
  targets=rows(out/'independent_h3_targets.jsonl');byid={r['id']:r for r in targets};need(len(byid)==len(targets)==8159,'abstract target coverage')
  flat(targets,out/'independent_h3_targets.flat')
  run('binary_sieve',[out/'independent_embedding_sieve',out/'independent_h3_targets.flat'],out/'independent_h3_binary_sieve.jsonl')
  py('verify_binary_witnesses.py');run('binary_absence_check',[out/'check_binary_absence'],out/'binary_absence_independent.json')
  need(json.loads((out/'binary_absence_independent.json').read_text())['hits']==0,'binary obstruction represented')
  binary=rows(out/'independent_h3_binary_sieve.jsonl');survivors=[byid[r['id']]for r in binary if r['status']=='survives']
  need(len(survivors)==1376,'binary survivor coverage')
  (out/'h3_source_inputs.gp').write_text('ids='+str([r['id']for r in survivors])+';\nHS=['+','.join(gpm(r['H'])for r in survivors)+'];\n')
  gpr('h3_source_orbits.gp');py('verify_source_orbits.py')
  reps=rows(out/'independent_h3_representatives.jsonl');rank={r['id']:r.get('short_rank')for r in binary}
  high=[r for r in reps if rank[r['id']]==4];low=[r for r in reps if rank[r['id']]<4]
  need(len(reps)==362 and len(high)==361 and len(low)==1,'source representative partition')
  highflat=out/'highspan.flat';flat(high,highflat)
  run('simultaneous_embedding',[out/'independent_simultaneous',highflat],out/'simultaneous.jsonl')
  py('verify_simultaneous.py','--targets',out/'independent_h3_targets.jsonl','--results',out/'simultaneous.jsonl','--quadruples',str(highflat)+'.baseline_quads','--output',out/'simultaneous_verification.json')
  py('test_affine_enum.py');py('prepare_lowspan_target.py')
  prepared=json.loads((out/'h3_lowspan_rebased.json').read_text());need(prepared['source_id']==low[0]['id'],'affine representative mismatch')
  py('affine_embedding_audit.py',out/'h3_lowspan_rebased.json','lowspan_affine')
  gpr('embedding_controls.gp');(out/'T64.flat').write_text((out/'embedding_controls.flat').read_text().splitlines()[0]+'\n')
  run('T64_positive_control',[out/'independent_simultaneous',out/'T64.flat'],out/'T64_positive.jsonl');need(rows(out/'T64_positive.jsonl')[0]['integral']>0,'positive embedding control failed')
  sim=rows(out/'simultaneous.jsonl');iv=json.loads((out/'simultaneous_verification.json').read_text());aff=json.loads((out/'lowspan_affine/summary.json').read_text());sm=json.loads((out/'h3_source_orbit_map.json').read_text())
  need(len(binary)==len({r['id']for r in binary})==len(targets)and {r['id']for r in binary}==set(byid),'binary full coverage')
  need(len(sm)==len({r[0]for r in sm})==len(survivors)and{r[0]for r in sm}=={r['id']for r in survivors},'source full coverage')
  need({r['id']for r in sim}=={r['id']for r in high}and len(sim)==len(high)and all(r['integral']==0 for r in sim),'simultaneous full coverage')
  need(iv['complete']and iv['source_frames']==len(high)and iv['quadruples_checked']==sum(r['quadruples']for r in sim)and iv['integral_embeddings']==0,'inverse verification incomplete')
  need(aff['completed']and aff['embeddings_found']==0 and[r['output_states']for r in aff['stages']]==[49,9,0]and all(r['independent_python']for r in aff['stages']),'affine exclusion incomplete')
  eliminated={r['id']for r in binary if r['status']!='survives'};er={r['id']for r in sim}|{low[0]['id']};eliminated|={t for t,rep,_ in sm if rep in er};need(eliminated==set(byid),'unaccounted target')
  usage=resource.getrusage(resource.RUSAGE_CHILDREN)
  summary.update(status='complete',rank4_h3_excluded=True,abstract_candidates=len(targets),binary_survivors=len(survivors),source_representatives=len(reps),highspan_representatives=len(high),lowspan_representatives=len(low),remaining=0,
   binary_status_counts={s:sum(r['status']==s for r in binary)for s in sorted({r['status']for r in binary})},simultaneous_counts={k:sum(r[k]for r in sim)for k in ['pairs','triples','quadruples','integral']},
   wall_seconds=time.perf_counter()-started,user_cpu_seconds=usage.ru_utime-usage0.ru_utime,system_cpu_seconds=usage.ru_stime-usage0.ru_stime)
  summary['output_hashes']={name:sha(out/name)for name in ['independent_h3_targets.jsonl','independent_h3_binary_sieve.jsonl','h3_source_orbit_map.json','simultaneous.jsonl','simultaneous_verification.json','h3_lowspan_rebased.json','lowspan_affine/summary.json']}
  save();print('H3_PROOF_COMPLETE',out/'FULL_RUN_SUMMARY.json',flush=True)
 except BaseException as e:
  summary.update(status='failed',error=str(e),wall_seconds=time.perf_counter()-started);save();raise
if __name__=='__main__':main()
