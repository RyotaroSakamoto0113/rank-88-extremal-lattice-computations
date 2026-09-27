#!/usr/bin/env python3
"""Verify one saved D4 shard without loading its production search functions."""
import argparse,csv,datetime,hashlib,json,re,resource,shutil,subprocess,sys,time
from pathlib import Path

def need(ok,s):
    if not ok:raise ValueError(s)

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def assignment_ids(path,prefix):
    text=path.read_text();ids=[int(x) for x in re.findall(r'^'+re.escape(prefix)+r'(\d+)\s*=',text,re.M)]
    need(len(ids)==len(set(ids)),f'duplicate {prefix} assignment');return ids

def gp_literal(x):
    if isinstance(x,(str,Path)):return json.dumps(str(x))
    if isinstance(x,list):return '['+','.join(gp_literal(v) for v in x)+']'
    return str(x)

def main():
    source=Path(__file__).resolve();root=source.parent.parent
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('shard',type=Path);p.add_argument('--output',type=Path);p.add_argument('--lookup',type=Path,default=root/'witnesses/ordinary_prune/ordinary_lookup.gp');p.add_argument('--exact-counts',type=Path);p.add_argument('--gp',default='gp');a=p.parse_args()
    shard=a.shard.resolve();out=a.output.resolve() if a.output else shard/'independent_certificate';out.mkdir(parents=True,exist_ok=False)
    counts_path=a.exact_counts.resolve() if a.exact_counts else shard/'exact_enumeration.stdout'
    required=['structure.gp','result.gp','exact_inputs.txt','exact_manifest.csv','enumerated_vectors.gp','quotient_metadata.gp','rank2_modules.gp','rank3_modules.gp','rank4_candidates.gp','pruning_witnesses.gp']
    for name in required:need((shard/name).is_file(),f'missing shard input {name}')
    need(counts_path.is_file(),'independent exact count file required')
    t0=time.perf_counter();c0=resource.getrusage(resource.RUSAGE_SELF);summary={'status':'running','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'shard':str(shard),'output':str(out),'input_hashes':{str(x):sha(x) for x in [*[shard/n for n in required],a.lookup,counts_path,source,source.with_suffix('.gp')]}}
    def save(): (out/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    save()
    try:
        rows=list(csv.reader((shard/'exact_manifest.csv').read_text().splitlines()));manifest=[]
        exact=[json.loads(x) for x in counts_path.read_text().splitlines() if x.strip()]
        need(len(rows)==len(exact)>0,'manifest and independent count coverage')
        for i,(row,ex) in enumerate(zip(rows,exact),1):
            stage={'second':2,'third':3,'fourth':4}.get(row[4]);need(stage is not None,'manifest stage');need(len(row)=={2:8,3:9,4:10}[stage],'manifest width')
            rr=[int(x) if j!=4 else stage for j,x in enumerate(row)];need(rr[0]==i,'contiguous manifest ID')
            need(int(ex['file'])==i and ex['dimension']==rr[1] and int(ex['bound'])==rr[2],'exact count form identity');need(ex['complete'] is True and ex['signed_vectors']==rr[3] and 2*ex['pairs']==rr[3],'complete independent exact counts match saved counts');manifest.append(rr)
        names={};files={'AUDIT_R_IDS':('rank2_modules.gp','R'),'AUDIT_T_IDS':('rank3_modules.gp','T'),'AUDIT_C_IDS':('rank4_candidates.gp','C'),'AUDIT_W2_IDS':('pruning_witnesses.gp','W2_'),'AUDIT_W3_IDS':('pruning_witnesses.gp','W3_')}
        for var,(filename,prefix) in files.items():names[var]=assignment_ids(shard/filename,prefix)
        for var in ['AUDIT_R_IDS','AUDIT_T_IDS','AUDIT_C_IDS']:need(names[var]==list(range(1,len(names[var])+1)),f'contiguous {var}')
        need(assignment_ids(shard/'quotient_metadata.gp','Q')==list(range(1,len(manifest)+1)),'quotient metadata covers all forms')
        need(assignment_ids(shard/'enumerated_vectors.gp','V')==list(range(1,len(manifest)+1)),'enumerated vectors cover all forms')
        with (shard/'exact_inputs.txt').open() as f, (out/'exact_inputs.gp').open('w') as dst:
            for row in manifest:
                line=next(f);n,b=map(int,line.split());need(n==row[1] and b==row[2],'integer input dimension and bound')
                A=[]
                for i in range(n):
                    ar=list(map(int,next(f).split()));need(len(ar)==n,'integer matrix row width');A.append(ar)
                dst.write('AUDIT_Q'+str(row[0])+'=['+';'.join(','.join(map(str,ar)) for ar in A)+'];\n')
            need(not any(line.strip() for line in f),'no extra exact forms')
        with (out/'metadata.gp').open('w') as f:
            f.write('AUDIT_MANIFEST='+gp_literal(manifest)+';\n')
            for var,ids in names.items():f.write(var+'='+gp_literal(ids)+';\n')
        lines=['default(nbthreads,1);','default(parisizemax,3000000000);','default(parisize,512000000);','AUDIT_OUT='+gp_literal(out)+';']
        for filename in ['structure.gp','result.gp','enumerated_vectors.gp','quotient_metadata.gp','rank2_modules.gp','rank3_modules.gp','rank4_candidates.gp','pruning_witnesses.gp']:lines.append('read('+gp_literal(shard/filename)+');')
        for path in [a.lookup,out/'metadata.gp',out/'exact_inputs.gp',source.with_suffix('.gp')]:lines.append('read('+gp_literal(path)+');')
        lines.append('quit;');(out/'launch.gp').write_text('\n'.join(lines)+'\n')
        started=datetime.datetime.now(datetime.timezone.utc).isoformat();t=time.perf_counter();c=resource.getrusage(resource.RUSAGE_CHILDREN)
        with (out/'verification.stdout').open('w') as stdout,(out/'verification.stderr').open('w') as stderr:
            result=subprocess.run([shutil.which(a.gp) or a.gp,'-q','-f',str(out/'launch.gp')],stdout=stdout,stderr=stderr)
        d=resource.getrusage(resource.RUSAGE_CHILDREN);summary['gp_stage']={'started_utc':started,'wall_seconds':time.perf_counter()-t,'user_cpu_seconds':d.ru_utime-c.ru_utime,'system_cpu_seconds':d.ru_stime-c.ru_stime,'returncode':result.returncode}
        need(result.returncode==0,'GP verifier exit code');errors=[x for x in (out/'verification.stderr').read_text().splitlines() if x.lstrip().startswith('***') and 'Warning:' not in x];need(not errors,'GP verifier errors: '+' / '.join(errors[:3]))
        text=(out/'verification.stdout').read_text();need('D4_INDEPENDENT_CERTIFICATE_VERIFIED' in text,'completion marker absent')
        rr=(out/'verification_result.gp').read_text();match=re.search(r'VERIFIED_COUNTS=(\[[^\]]*\]);VERIFIED_ANCHOR=(\d+);',rr);need(match is not None,'verified results absent')
        values=json.loads(match[1]);summary.update(status='verified',anchor=int(match[2]),counts=values,independent_exact_forms=len(exact),all_raw_vectors_covered_once=True,all_pruning_witnesses_valid=True,all_quotient_forms_reconstructed=True,independent_complete_counts_match=True,no_rank4_candidates=(values[11]==0))
    except BaseException as e:
        summary.update(status='failed',error=type(e).__name__+': '+str(e));raise
    finally:
        c1=resource.getrusage(resource.RUSAGE_SELF);summary['total_wall_seconds']=time.perf_counter()-t0;summary['python_user_cpu_seconds']=c1.ru_utime-c0.ru_utime;summary['python_system_cpu_seconds']=c1.ru_stime-c0.ru_stime;summary['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();save()
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
