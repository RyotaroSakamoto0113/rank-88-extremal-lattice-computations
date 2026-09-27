#!/usr/bin/env python3
"""Generate and independently count all quotient spheres for projective d3<=26."""
import datetime,hashlib,json,pathlib,re,resource,shutil,subprocess,time,tempfile
HERE=pathlib.Path(__file__).resolve().parent
def need(b,s):
    if not b:raise RuntimeError(s)
def main():
    boot=pathlib.Path((HERE/'LATEST_BOOTSTRAP.txt').read_text().strip());need(json.loads((boot/'BOOTSTRAP_SUMMARY.json').read_text())['status']=='complete','bootstrap incomplete')
    root=pathlib.Path(tempfile.mkdtemp(prefix=datetime.datetime.now().strftime('%Y%m%d-%H%M%S-d3-'),dir=HERE/'runs'));out=root/'d3';out.mkdir();(root/'src').mkdir()
    for f in(HERE/'src').glob('*.gp'):shutil.copy2(f,root/'src'/f.name)
    for f in ['ambient.gp','generators.gp']:shutil.copy2(boot/'bootstrap'/f,out/f)
    launcher=root/'d3.gp';launcher.write_text('default(parisizemax,3000000000);\ndefault(nbthreads,1);\nOUT='+json.dumps(str(out))+';BOOT='+json.dumps(str(boot/'bootstrap'))+';\nread('+json.dumps(str(root/'src/common.gp'))+');\nread('+json.dumps(str(root/'src/d3.gp'))+');\n')
    summary={'status':'running','run_directory':str(root),'bootstrap_directory':str(boot),'stages':[]};start=time.perf_counter();u0=resource.getrusage(resource.RUSAGE_CHILDREN)
    def save():(root/'D3_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    def run(name,cmd):
        print('START',name,flush=True);t=time.perf_counter();u=resource.getrusage(resource.RUSAGE_CHILDREN)
        with(root/(name+'.stdout')).open('w')as a,(root/(name+'.stderr')).open('w')as b:r=subprocess.run(list(map(str,cmd)),stdout=a,stderr=b)
        v=resource.getrusage(resource.RUSAGE_CHILDREN);summary['stages'].append({'stage':name,'wall_seconds':time.perf_counter()-t,'user_cpu_seconds':v.ru_utime-u.ru_utime,'system_cpu_seconds':v.ru_stime-u.ru_stime,'returncode':r.returncode});save();need(r.returncode==0,name+' failed')
        need(not any(x.lstrip().startswith('***')and'Warning:'not in x for x in(root/(name+'.stderr')).read_text().splitlines()),name+' PARI error');print('DONE',name,summary['stages'][-1]['wall_seconds'],flush=True)
    print('RUN_DIRECTORY',root,flush=True)
    try:
        run('quotient_generation',['gp','-fq',launcher]);need('D3_ALL_PROJECTIVE_COMPLETE'in(root/'quotient_generation.stdout').read_text(),'generation incomplete')
        run('exact_quotient_counts',[boot/'exact_enum','--batch',out/'exact_inputs.txt','60'])
        rr=[json.loads(x)for x in(root/'exact_quotient_counts.stdout').read_text().splitlines()if x.strip()];mm=[x.split(',')for x in(out/'exact_manifest.csv').read_text().splitlines()if x.strip()]
        need(len(rr)==len(mm)>0,'batch count');need(all(r['complete']and r['signed_vectors']==int(m[3])for r,m in zip(rr,mm)),'exact count mismatch/incomplete')
        result=json.loads(re.search(r'RESULT=(\[.*?\]);',(out/'result.gp').read_text()).group(1));need(result[-1]==len(mm),'GP/batch coverage')
        labels=['first_orbits','raw_rank2','nonsaturated_rank2','shorter_line_rejections','duplicate_rank2','kept_rank2','raw_rank3','nonsaturated_rank3','saturated_rank3','free_rank3','min_rank2','min_rank3','exact_forms']
        summary['counts']=dict(zip(labels,result));summary['d1_projective_lower']=6;summary['d2_projective_lower']=16;summary['d3_projective_lower']=24;summary['d3_free_lower']=27
        need(result[9]==0 and result[10]>=16 and result[11]>=24,'claimed lower bound contradicted')
        u=resource.getrusage(resource.RUSAGE_CHILDREN);summary.update(status='complete',wall_seconds=time.perf_counter()-start,user_cpu_seconds=u.ru_utime-u0.ru_utime,system_cpu_seconds=u.ru_stime-u0.ru_stime)
        summary['input_hashes']={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest()for f in(root/'src').glob('*')};save();(HERE/'LATEST_D3.txt').write_text(str(root)+'\n');print('D1_D2_D3_LOWER_BOUNDS_CERTIFIED',summary['wall_seconds'],flush=True)
    except BaseException as e:summary.update(status='failed',error=str(e),wall_seconds=time.perf_counter()-start);save();raise
if __name__=='__main__':main()
