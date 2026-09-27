"""Regenerate the h=3 candidates and sieve against J0, retaining all evidence."""
import collections,json,pathlib,re,subprocess,sys,time
HERE=pathlib.Path(__file__).resolve().parent;ROOT=HERE/'run';OUT=ROOT/'output/rank4_audit'
(OUT/'logs').mkdir(exist_ok=True)
timings=[]
def run(label,cmd,dest=None):
 print('START',label,flush=True);t=time.perf_counter()
 f=pathlib.Path(dest)if dest else OUT/'logs'/(label+'.stdout');err=OUT/'logs'/(label+'.stderr')
 with f.open('w')as fo,err.open('w')as fe:r=subprocess.run(list(map(str,cmd)),cwd=ROOT,stdout=fo,stderr=fe)
 timings.append({'stage':label,'seconds':time.perf_counter()-t,'returncode':r.returncode})
 (OUT/'initial_timings.json').write_text(json.dumps(timings,indent=2)+'\n')
 assert r.returncode==0,(label,err.read_text())
 if 'gp'==pathlib.Path(cmd[0]).name:assert not any(l.lstrip().startswith('***')and'Warning:'not in l for l in err.read_text().splitlines()),err.read_text()
 print('DONE',label,round(timings[-1]['seconds'],3),flush=True)
def rows(p):return [json.loads(l)for l in p.read_text().splitlines()if l.strip()]
def flat(ts,p):p.write_text(''.join(str(t['id'])+' '+' '.join(str(v)for row in t['H']for z in row for v in z)+'\n'for t in ts))
def gpm(H):return '['+';'.join(','.join(f'({a})+({b})*w'for a,b in row)for row in H)+']'
def main():
 for n,gmp in [('abstract_h3',False),('check_M_shell_orbits',False),('check_binary_absence',False),('independent_embedding_sieve',True),('independent_simultaneous',True),('exact_affine_enum',True)]:
  run('compile_'+n,['clang++','-O3','-std=c++17',OUT/(n+'.cpp'),'-o',OUT/n,*(['-I/opt/homebrew/include','-L/opt/homebrew/lib','-lgmpxx','-lgmp']if gmp else[])])
 exact=ROOT/'output/d3/benchmark/exact_enum'
 run('compile_exact_enum',['clang++','-O3','-std=c++17',str(exact)+'.cpp','-o',exact,'-I/opt/homebrew/include','-L/opt/homebrew/lib','-lgmpxx','-lgmp'])
 run('J0_ambient',['gp','-fq',HERE/'ambient.gp'])
 for n in ['M_integer_data.txt','M_shell12.txt','M_shell12_count_input.txt']:
  p=OUT/n;p.write_text(re.sub(r'[\[\],;]',' ',p.read_text()))
 run('J0_exact_sphere_count',[exact,OUT/'M_shell12_count_input.txt','60'],OUT/'M_shell12_exact_count.json')
 run('J0_shell_orbits',[OUT/'check_M_shell_orbits'],OUT/'M_shell12_orbit_audit.json')
 ec=json.loads((OUT/'M_shell12_exact_count.json').read_text());orb=json.loads((OUT/'M_shell12_orbit_audit.json').read_text())
 assert ec['complete']and ec['signed_vectors']==orb['signed_count']==2*orb['pair_count']
 assert orb['unique_valid']and orb['generator_closed']and all(s['q']>=6 and s['q']!=7 for s in orb['shells'])
 print('SHELL',json.dumps(orb),flush=True)
 d=(OUT/'M_integer_data.txt').read_text().split();(OUT/'h2_M_input.txt').write_text('22 18\n'+' '.join(d[1:1+2*22*22])+'\n')
 run('charts',['gp','-fq',OUT/'independent_chart_audit.gp']);run('export_charts',['gp','-fq',OUT/'export_charts_json.gp'])
 run('abstract',[OUT/'abstract_h3',OUT/'independent_h3','27','48','1'])
 ts=rows(OUT/'independent_h3_targets.jsonl');assert len(ts)==8159
 flat(ts,OUT/'independent_h3_targets.flat')
 run('binary_sieve',[OUT/'independent_embedding_sieve',OUT/'independent_h3_targets.flat'],OUT/'independent_h3_binary_sieve.jsonl')
 run('binary_witness_verification',[sys.executable,OUT/'verify_binary_witnesses.py'])
 run('binary_absence_verification',[OUT/'check_binary_absence'],OUT/'binary_absence_independent.json')
 bs=rows(OUT/'independent_h3_binary_sieve.jsonl');ids={r['id']for r in bs if r['status']=='survives'};ss=[t for t in ts if t['id']in ids]
 print('BINARY',dict(collections.Counter(r['status']for r in bs)),flush=True)
 (OUT/'h3_source_inputs.gp').write_text('ids='+str([r['id']for r in ss])+';\nHS=['+','.join(gpm(r['H'])for r in ss)+'];\n')
 run('source_orbits',['gp','-fq',OUT/'h3_source_orbits.gp']);run('verify_source_orbits',[sys.executable,OUT/'verify_source_orbits.py'])
 reps=rows(OUT/'independent_h3_representatives.jsonl');rank={r['id']:r.get('short_rank')for r in bs}
 high=[r for r in reps if rank[r['id']]==4];low=[r for r in reps if rank[r['id']]<4]
 flat(high,OUT/'highspan.flat');(OUT/'lowspan_targets.json').write_text(json.dumps(low,indent=2)+'\n')
 print('REPRESENTATIVES',len(reps),'HIGH',len(high),'LOW',[(r['id'],rank[r['id']])for r in low],flush=True)
 run('simultaneous',[OUT/'independent_simultaneous',OUT/'highspan.flat'],OUT/'simultaneous.jsonl')
 run('verify_simultaneous',[sys.executable,OUT/'verify_simultaneous.py','--targets',OUT/'independent_h3_targets.jsonl','--results',OUT/'simultaneous.jsonl','--quadruples',str(OUT/'highspan.flat')+'.baseline_quads','--output',OUT/'simultaneous_verification.json'])
 sim=rows(OUT/'simultaneous.jsonl');print('SIMULTANEOUS', {k:sum(r[k]for r in sim)for k in ['pairs','triples','quadruples','integral']},flush=True)
 result={'abstract_candidates':len(ts),'binary_counts':dict(collections.Counter(r['status']for r in bs)),'representatives':len(reps),'highspan':len(high),'lowspan':len(low),'integral_embeddings':sum(r['integral']for r in sim)}
 (OUT/'INITIAL_RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
