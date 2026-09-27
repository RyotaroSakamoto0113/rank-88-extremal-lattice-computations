#!/usr/bin/env python3
"""Fresh q<=12 sphere, exact count and complete ordinary orbit reconstruction."""
import datetime,hashlib,json,pathlib,re,resource,shutil,subprocess,time,tempfile
HERE=pathlib.Path(__file__).resolve().parent;PROJECT=HERE.parents[1]
def need(b,s):
    if not b:raise RuntimeError(s)
def main():
    boot=pathlib.Path((HERE/'LATEST_BOOTSTRAP.txt').read_text().strip());root=pathlib.Path(tempfile.mkdtemp(prefix=datetime.datetime.now().strftime('%Y%m%d-%H%M%S-ordinary-'),dir=HERE/'runs'));out=root/'output/rank4_audit';out.mkdir(parents=True);(root/'output/d3/benchmark').mkdir(parents=True)
    shutil.copy2(PROJECT/'output/d3/benchmark/ambient.gp',root/'output/d3/benchmark/ambient.gp');shutil.copy2(PROJECT/'output/rank4_audit/ambient_qfauto.gp',out/'ambient_qfauto.gp');shutil.copy2(PROJECT/'output/rank4_d27_d48/src/independent_M_shell.gp',out/'independent_M_shell.gp')
    t=time.perf_counter();u0=resource.getrusage(resource.RUSAGE_CHILDREN);s={'status':'running','run_directory':str(root),'stages':[]}
    def save():(root/'ORDINARY_SUMMARY.json').write_text(json.dumps(s,indent=2)+'\n')
    def run(name,cmd):
        print('START',name,flush=True);tt=time.perf_counter();u=resource.getrusage(resource.RUSAGE_CHILDREN)
        with(root/(name+'.stdout')).open('w')as a,(root/(name+'.stderr')).open('w')as b:r=subprocess.run(list(map(str,cmd)),cwd=root,stdout=a,stderr=b)
        v=resource.getrusage(resource.RUSAGE_CHILDREN);s['stages'].append({'stage':name,'wall_seconds':time.perf_counter()-tt,'user_cpu_seconds':v.ru_utime-u.ru_utime,'system_cpu_seconds':v.ru_stime-u.ru_stime,'returncode':r.returncode});save();need(r.returncode==0,name+' failed');need(not any(x.lstrip().startswith('***')and'Warning:'not in x for x in(root/(name+'.stderr')).read_text().splitlines()),name+' PARI error');print('DONE',name,s['stages'][-1]['wall_seconds'],flush=True)
    print('RUN_DIRECTORY',root,flush=True)
    try:
        run('sphere_generation',['gp','-fq',out/'independent_M_shell.gp'])
        for name in ['M_integer_data.txt','M_shell12.txt','M_shell12_count_input.txt']:
            f=out/name;f.write_text(re.sub(r'[\[\],;]',' ',f.read_text()))
        run('exact_sphere_count',[boot/'exact_enum',out/'M_shell12_count_input.txt','60']);e=json.loads((root/'exact_sphere_count.stdout').read_text());need(e['complete']and e['signed_vectors']==567226,'sphere incomplete')
        run('orbits_stabilizers',[boot/'minima_orbits',out/'M_integer_data.txt',out/'M_shell12.txt',out/'orbits','24']);o=json.loads((root/'orbits_stabilizers.stdout').read_text());need(o['status']=='verified'and o['group_order']==12144 and o['orbit_count']==78 and o['pair_count']*2==e['signed_vectors'],'orbits incomplete')
        for name in ['point_orbit_lookup.txt','orbit_representatives.gp','stabilizers.jsonl','group_matrices.txt']:
            need((out/'orbits'/name).read_bytes()==(HERE/'witnesses/orbits'/name).read_bytes(),'fresh orbit data differs '+name)
        u=resource.getrusage(resource.RUSAGE_CHILDREN);s.update(status='complete',wall_seconds=time.perf_counter()-t,user_cpu_seconds=u.ru_utime-u0.ru_utime,system_cpu_seconds=u.ru_stime-u0.ru_stime,signed_vectors=e['signed_vectors'],ordinary_orbits=78,previous_lookup_data_identical=True);save();(HERE/'LATEST_ORDINARY.txt').write_text(str(root)+'\n');print('ORDINARY_SPHERE_CERTIFIED',s['wall_seconds'],flush=True)
    except BaseException as e:s.update(status='failed',error=str(e),wall_seconds=time.perf_counter()-t);save();raise
if __name__=='__main__':main()
