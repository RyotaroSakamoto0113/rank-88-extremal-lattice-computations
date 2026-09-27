#!/usr/bin/env python3
"""Link the actual D4 snapshot to certified first modules, stabilizers and lookups."""
from pathlib import Path
import hashlib,json,subprocess,time,resource
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
RUN=Path((ROOT/'LATEST_D4_RUN.txt').read_text().strip());OUT=HERE/'d4_snapshot_links';OUT.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(b,s):
 if not b:raise RuntimeError(s)
initial=json.loads((ROOT/'initial4/verification.json').read_text());stabs=json.loads((ROOT/'initial4/stabilizers_verification.json').read_text())
inputs=RUN/'inputs';links=[]
for name in ['FIRST4.gp','STAB4.gp']:
 src=ROOT/'initial4'/name;dst=inputs/name;expected=initial['sha256'][name]if name=='FIRST4.gp'else stabs['sha256'][str(src)]
 need(sha(src)==sha(dst)==expected,'snapshot certified hash mismatch '+name);links.append({'name':name,'sha256':expected})
for src,dst in [(ROOT/'witnesses/orbits/orbit_representatives.gp',inputs/'ordinary_representatives.gp'),(ROOT/'witnesses/ordinary_prune/ordinary_lookup.gp',inputs/'ordinary_prune/ordinary_lookup.gp')]:
 need(sha(src)==sha(dst),'snapshot source mismatch '+src.name);links.append({'name':src.name,'sha256':sha(src)})
launch=OUT/'compare.gp';snapfirst=OUT/'runtime_first.json'
launch.write_text('default(parisizemax,3000000000);\ndefault(parisize,100000000);\ndefault(nbthreads,1);\nread('+json.dumps(str(inputs/'FIRST4.gp'))+');\nread('+json.dumps(str(inputs/'STAB4.gp'))+');\nread('+json.dumps(str(inputs/'ordinary_representatives.gp'))+');\nneed(b,s)={if(!b,error(s))};\nrows(A)=vector(matsize(A)[1],i,vector(matsize(A)[2],j,A[i,j]));\nwrite('+json.dumps(str(snapfirst))+',vector(#FIRST4,i,[FIRST4[i][1],FIRST4[i][2],FIRST4[i][3],FIRST4[i][4],rows(FIRST4[i][5])]));\nfor(i=1,#STAB4,x=MINIMA_ORBIT_REPS[i][4]~;ss=Set(STAB4[i]);need(#ss==#STAB4[i],"distinct stabilizer matrices");need(setsearch(ss,matid(22)),"stabilizer identity");for(j=1,#ss,need(ss[j]*x==x||ss[j]*x==-x,"actual signed vector stabilized");for(k=1,#ss,need(setsearch(ss,ss[j]*ss[k]),"stabilizer group closure"))));\nprint("SNAPSHOT_STABILIZERS_CLOSED_AND_FIX_SIGNED_X");quit;\n')
snapfirst.unlink(missing_ok=True);u=resource.getrusage(resource.RUSAGE_CHILDREN);t=time.perf_counter();p=subprocess.run(['gp','-fq',str(launch)],capture_output=True,text=True);wall=time.perf_counter()-t;v=resource.getrusage(resource.RUSAGE_CHILDREN)
(OUT/'compare.stdout').write_text(p.stdout);(OUT/'compare.stderr').write_text(p.stderr)
need(p.returncode==0 and 'SNAPSHOT_STABILIZERS_CLOSED_AND_FIX_SIGNED_X'in p.stdout and not any(x.lstrip().startswith('***')and'Warning:'not in x for x in p.stderr.splitlines()),'GP snapshot verification failed')
actual=json.loads(snapfirst.read_text());records=[json.loads(x)for x in(ROOT/'initial4/records.jsonl').read_text().splitlines()if x.strip()]
need(actual==[[r[0],r[1],r[2],r[3],r[8]]for r in records],'FIRST4 does not encode independently verified records')
summary={'status':'verified','run_directory':str(RUN),'first_module_count':len(actual),'semantic_initial_records_match':True,'each_stabilizer_set_is_closed_group':True,'each_stabilizer_fixes_actual_signed_anchor':True,'snapshot_links':links,'gp_wall_seconds':wall,'gp_user_cpu_seconds':v.ru_utime-u.ru_utime,'gp_system_cpu_seconds':v.ru_stime-u.ru_stime}
(OUT/'verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
