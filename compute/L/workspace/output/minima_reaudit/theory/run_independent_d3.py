#!/usr/bin/env python3
import pathlib,json,time,resource,subprocess,sys,hashlib,tempfile,shutil
HERE=pathlib.Path(__file__).resolve().parent;MAIN=HERE.parent
boot=pathlib.Path((MAIN/'LATEST_BOOTSTRAP.txt').read_text().strip());root=pathlib.Path(tempfile.mkdtemp(prefix='independent-d3-',dir=MAIN/'runs'));out=root/'data';out.mkdir()
shutil.copy2(HERE/'independent_d3.gp',root/'independent_d3.gp')
launcher=root/'launch.gp';launcher.write_text('default(parisizemax,3000000000);\nBOOT='+json.dumps(str(boot/'bootstrap'))+';OUT='+json.dumps(str(out))+';\nread('+json.dumps(str(root/'independent_d3.gp'))+');\n')
summary={'status':'running','bootstrap_directory':str(boot),'run_directory':str(root),'stages':[],'algorithm':'Independent qflllgram kernel/image quotient; all binary modules; no short-line pruning'};start=time.perf_counter()
def save():(root/'INDEPENDENT_D3_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
def need(b,s):
 if not b:raise RuntimeError(s)
def run(name,cmd):
 t=time.perf_counter();u=resource.getrusage(resource.RUSAGE_CHILDREN)
 with (root/(name+'.stdout')).open('w')as so,(root/(name+'.stderr')).open('w')as se:r=subprocess.run(list(map(str,cmd)),stdout=so,stderr=se)
 v=resource.getrusage(resource.RUSAGE_CHILDREN);summary['stages'].append({'name':name,'wall_seconds':time.perf_counter()-t,'user_cpu_seconds':v.ru_utime-u.ru_utime,'system_cpu_seconds':v.ru_stime-u.ru_stime,'returncode':r.returncode});save();need(r.returncode==0,name+' failed');need(not any(s.lstrip().startswith('***')and'Warning:'not in s for s in(root/(name+'.stderr')).read_text().splitlines()),name+' GP error');print(name,summary['stages'][-1],flush=True)
print('RUN_DIRECTORY',root,flush=True);save()
try:
 run('generate',['gp','-fq',launcher]);need('INDEPENDENT_D3_ENUMERATION_COMPLETE'in(root/'generate.stdout').read_text(),'missing completion')
 run('exact_counts',[boot/'exact_enum','--batch',out/'exact_inputs.txt','60'])
 rr=[json.loads(x)for x in(root/'exact_counts.stdout').read_text().splitlines()if x.strip()];mm=[x.split(',')for x in(out/'exact_manifest.csv').read_text().splitlines()if x.strip()]
 need(len(rr)==len(mm)>0,'coverage');need(all(r['complete']and r['signed_vectors']==int(m[3])for r,m in zip(rr,mm)),'exact count mismatch/incomplete')
 summary['counts']=json.loads((out/'independent_result.json').read_text());need(summary['counts']['forms']==len(mm),'batch count');need(summary['counts']['min_binary']>=16 and summary['counts']['free_ternary']==0 and summary['counts']['min_ternary']>=24,'lower bounds contradicted')
 summary.update(status='complete',wall_seconds=time.perf_counter()-start,input_hashes={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest()for f in[root/'independent_d3.gp',launcher]},all_exact_forms=len(rr));save();(HERE/'LATEST_INDEPENDENT_D3.txt').write_text(str(root)+'\n');print('INDEPENDENT_D3_CERTIFIED',summary,flush=True)
except BaseException as e:summary.update(status='failed',wall_seconds=time.perf_counter()-start,error=str(e));save();raise
