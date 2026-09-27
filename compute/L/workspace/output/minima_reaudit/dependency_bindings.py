#!/usr/bin/env python3
"""Bind every successful computation to the same audited J and initial data."""
import datetime, hashlib, json, pathlib, resource, subprocess, time
HERE = pathlib.Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def need(x, s):
    if not x: raise ValueError(s)
def latest(name):
    return pathlib.Path((HERE/name).read_text().strip())
def main():
    start=time.perf_counter(); cpu=resource.getrusage(resource.RUSAGE_SELF)
    boot=latest('LATEST_BOOTSTRAP.txt'); d3=latest('LATEST_D3.txt'); ordinary=latest('LATEST_ORDINARY.txt'); d4=latest('LATEST_D4.txt')
    indep=pathlib.Path((HERE/'theory/LATEST_INDEPENDENT_D3.txt').read_text().strip())
    checks=[]
    def same(a,b):
        need(sha(a)==sha(b),'different dependency: '+str(a)+' / '+str(b))
        checks.append({'left':str(a),'right':str(b),'sha256':sha(a)})
    for p,n in [(boot,'BOOTSTRAP_SUMMARY.json'),(d3,'D3_SUMMARY.json'),(ordinary,'ORDINARY_SUMMARY.json'),(d4,'D4_SUMMARY.json')]:
        need(read(p/n)['status']=='complete','incomplete '+str(p))
    need(read(indep/'DERIVED_LOWER_BOUNDS.json')['d3_projective_lower']==27,'projective d3 prerequisite')
    need('AMBIENT_INPUT_VERIFIED' in (boot/'input_identity.stdout').read_text(),'original J identity')
    manifest=read(d4/'INPUTS_SHA256.json')
    for name,h in manifest.items(): need(sha(d4/name)==h,'changed run snapshot '+name)
    for name in ['ambient.gp','generators.gp']:
        same(d4/'inputs'/name,boot/'bootstrap'/name)
        same(d4/'inputs'/name,d3/'d3'/name)
    same(boot/'bootstrap/ambient.gp',PROJECT/'output/d3/benchmark/ambient.gp')
    same(d4/'inputs/FIRST4.gp',HERE/'initial4/FIRST4.gp')
    need(sha(d4/'inputs/FIRST4.gp')==read(HERE/'initial4/verification.json')['sha256']['FIRST4.gp'],'initial class certificate hash')
    same(d4/'inputs/STAB4.gp',HERE/'initial4/STAB4.gp')
    stab=read(HERE/'initial4/stabilizers_verification.json')
    need(any(pathlib.Path(k).name=='STAB4.gp' and v==sha(d4/'inputs/STAB4.gp') for k,v in stab['sha256'].items()),'stabilizer certificate hash')
    orb=ordinary/'output/rank4_audit/orbits'
    for n in ['point_orbit_lookup.txt','orbit_representatives.gp','stabilizers.jsonl','group_matrices.txt']:
        same(orb/n,HERE/'witnesses/orbits'/n)
    same(d4/'inputs/ordinary_representatives.gp',orb/'orbit_representatives.gp')
    same(d4/'inputs/ordinary_prune/ordinary_lookup.gp',HERE/'witnesses/ordinary_prune/ordinary_lookup.gp')
    lookup=read(HERE/'witnesses/ordinary_prune/lookup_conversion.json')
    need(lookup['output']['sha256']==sha(d4/'inputs/ordinary_prune/ordinary_lookup.gp'),'lookup conversion output')
    need(any(pathlib.Path(k).name=='point_orbit_lookup.txt' and v==sha(orb/'point_orbit_lookup.txt') for k,v in lookup['inputs'].items()),'lookup conversion source')
    need(list(map(int,(ordinary/'output/rank4_audit/M_integer_data.txt').read_text().split()))==list(map(int,(HERE/'witnesses/M_integer_data.txt').read_text().split())),'upper and lower input matrices')
    same(d4/'src/common.gp',d3/'src/common.gp')
    same(d4/'src/common.gp',boot/'src/common.gp')
    out=HERE/'dependency_checks';out.mkdir(exist_ok=True)
    def gpstr(p): return json.dumps(str(p))
    lines=['default(nbthreads,1);default(parisizemax,3000000000);',
           'read('+gpstr(d4/'inputs/ambient.gp')+');GREF=G;WREF=W;',
           'read('+gpstr(boot/'bootstrap/cusps.gp')+');CCREF=CC;',
           'read('+gpstr(d4/'inputs/FIRST4.gp')+');',
           'read('+gpstr(d4/'inputs/STAB4.gp')+');',
           'read('+gpstr(d4/'inputs/ordinary_representatives.gp')+');',
           'need(b,s)=if(!b,error(s));']
    for i in range(78):
        shard=d4/'anchors'/f'{i:03d}'
        need(read(shard/'ANCHOR_SUMMARY.json')['status']=='complete','incomplete shard')
        same(shard/'ambient.gp',d4/'inputs/ambient.gp')
        same(shard/'generators.gp',d4/'inputs/generators.gp')
        lines += ['read('+gpstr(shard/'structure.gp')+');',
          f'need(SAVED_G==GREF && SAVED_W==WREF && SAVED_CC==CCREF && SAVED_NN==[1,2,2],"shard {i} ambient/cusps");',
          f'need(SAVED_FIRST==FIRST4[{i+1}] && SAVED_STAB==STAB4[{i+1}] && SAVED_X==MINIMA_ORBIT_REPS[{i+1}][4]~,"shard {i} initial/stabilizer");']
    lines += ['print("ALL_78_DEPENDENCY_BINDINGS_VERIFIED");quit;']
    (out/'check.gp').write_text('\n'.join(lines)+'\n')
    t=time.perf_counter(); c=resource.getrusage(resource.RUSAGE_CHILDREN)
    with (out/'check.stdout').open('w') as a, (out/'check.stderr').open('w') as b:
        p=subprocess.run(['gp','-fq',str(out/'check.gp')],stdout=a,stderr=b)
    v=resource.getrusage(resource.RUSAGE_CHILDREN)
    need(p.returncode==0 and 'ALL_78_DEPENDENCY_BINDINGS_VERIFIED' in (out/'check.stdout').read_text(),'GP binding verification failed')
    need(not any(x.lstrip().startswith('***') and 'Warning:' not in x for x in (out/'check.stderr').read_text().splitlines()),'GP binding error')
    result={'status':'verified','d4_run':str(d4),'initial_classes_and_stabilizers_match_verified_inputs':True,'all_78_shard_structures_match':True,'upper_and_lower_bounds_use_same_J':True,'checks':checks,'gp_stage':{'wall_seconds':time.perf_counter()-t,'user_cpu_seconds':v.ru_utime-c.ru_utime,'system_cpu_seconds':v.ru_stime-c.ru_stime},'wall_seconds':time.perf_counter()-start,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    result['python_cpu_seconds']=resource.getrusage(resource.RUSAGE_SELF).ru_utime-cpu.ru_utime
    (HERE/'dependency_bindings.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
if __name__=='__main__':main()
