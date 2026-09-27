#!/usr/bin/env python3
"""Reconstruct and certify the three small projective-line spheres."""
import datetime,hashlib,json,os,pathlib,resource,shutil,subprocess,sys,tempfile,time
HERE=pathlib.Path(__file__).resolve().parent;PROJECT=HERE.parents[1]
def need(b,s):
    if not b:raise RuntimeError(s)
def main():
    root=pathlib.Path(tempfile.mkdtemp(prefix=datetime.datetime.now().strftime('%Y%m%d-%H%M%S-'),dir=HERE/'runs'));out=root/'bootstrap';out.mkdir();(root/'src').mkdir()
    for f in (HERE/'src').glob('*.gp'):shutil.copy2(f,root/'src'/f.name)
    for source,name in [('output/d3/benchmark/ambient.gp','ambient.gp'),('output/rank4_audit/ambient_qfauto.gp','generators.gp')]:shutil.copy2(PROJECT/source,out/name)
    shutil.copy2(HERE/'witnesses/minima_orbits.cpp',root/'src/minima_orbits.cpp');shutil.copy2(PROJECT/'output/d3/benchmark/exact_enum.cpp',root/'src/exact_enum.cpp')
    summary={'status':'running','run_directory':str(root),'stages':[]};t=time.perf_counter();u0=resource.getrusage(resource.RUSAGE_CHILDREN)
    def save(): (root/'BOOTSTRAP_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    def run(name,cmd,cwd=PROJECT):
        print('START',name,flush=True);u=resource.getrusage(resource.RUSAGE_CHILDREN);st=time.perf_counter()
        with (root/(name+'.stdout')).open('w')as a,(root/(name+'.stderr')).open('w')as b:r=subprocess.run(list(map(str,cmd)),cwd=cwd,stdout=a,stderr=b)
        v=resource.getrusage(resource.RUSAGE_CHILDREN);summary['stages'].append({'stage':name,'wall_seconds':time.perf_counter()-st,'user_cpu_seconds':v.ru_utime-u.ru_utime,'system_cpu_seconds':v.ru_stime-u.ru_stime,'returncode':r.returncode});save()
        need(r.returncode==0,name+' failed');need(not any(x.lstrip().startswith('***')and'Warning:'not in x for x in(root/(name+'.stderr')).read_text().splitlines()),name+' PARI error')
        print('DONE',name,summary['stages'][-1]['wall_seconds'],flush=True)
    print('RUN_DIRECTORY',root,flush=True)
    try:
        run('input_identity',['gp','-fq',PROJECT/'output/rank4_audit/verify_ambient_input.gp'])
        need('AMBIENT_INPUT_VERIFIED'in(root/'input_identity.stdout').read_text(),'input identity incomplete')
        launcher=root/'bootstrap.gp';launcher.write_text('default(parisizemax,3000000000);\ndefault(nbthreads,1);\nOUT='+json.dumps(str(out))+';\nread('+json.dumps(str(root/'src/common.gp'))+');\nread('+json.dumps(str(root/'src/bootstrap.gp'))+');\n')
        run('bootstrap_gp',['gp','-fq',launcher]);need('BOOTSTRAP_COMPLETE'in(root/'bootstrap_gp.stdout').read_text(),'bootstrap not complete')
        run('compile_orbits',['clang++','-O3','-std=c++17',root/'src/minima_orbits.cpp','-o',root/'minima_orbits'])
        run('compile_exact',['clang++','-O3','-std=c++17','-I/opt/homebrew/include','-L/opt/homebrew/lib',root/'src/exact_enum.cpp','-lgmpxx','-lgmp','-o',root/'exact_enum'])
        for c in range(1,4):run('orbits_'+str(c),[root/'minima_orbits',out/f'cusp{c}_integer_data.txt',out/f'cusp{c}_shell.txt',out/f'orbits{c}','20','--orbits-only'])
        run('exact_spheres',[root/'exact_enum','--batch',out/'exact_inputs.txt','60'])
        rr=[json.loads(x)for x in(root/'exact_spheres.stdout').read_text().splitlines()if x.strip()];mm=[x.split(',')for x in(out/'exact_manifest.csv').read_text().splitlines()if x.strip()]
        need(len(rr)==len(mm)==3,'sphere batch coverage')
        for r,m in zip(rr,mm):need(r['complete']and r['signed_vectors']==int(m[3]),'sphere incomplete/count mismatch')
        for c in range(1,4):
            a=json.loads((root/f'orbits_{c}.stdout').read_text());need(a['status']=='verified','orbit computation incomplete');summary.setdefault('cusps',[]).append(a)
        u=resource.getrusage(resource.RUSAGE_CHILDREN);summary.update(status='complete',wall_seconds=time.perf_counter()-t,user_cpu_seconds=u.ru_utime-u0.ru_utime,system_cpu_seconds=u.ru_stime-u0.ru_stime)
        summary['input_hashes']={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest()for f in(root/'src').glob('*')};save();(HERE/'LATEST_BOOTSTRAP.txt').write_text(str(root)+'\n');print('BOOTSTRAP_CERTIFIED',summary['wall_seconds'],flush=True)
    except BaseException as e:summary.update(status='failed',error=str(e),wall_seconds=time.perf_counter()-t);save();raise
if __name__=='__main__':main()
