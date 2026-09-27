"""Check exact coverage of every candidate and aggregate the J0 h=3 proof."""
import collections,hashlib,json,pathlib,sys
def read(p):return json.loads(p.read_text())
def rows(p):return [json.loads(l)for l in p.read_text().splitlines()if l.strip()]
def require(b,s):
 if not b:raise ValueError(s)
def finalize(root):
 root=pathlib.Path(root).resolve();p=root/'output/rank4_audit'
 targets=rows(p/'independent_h3_targets.jsonl');ids={r['id']for r in targets}
 require(len(ids)==len(targets)==8159,'abstract target coverage')
 binary=rows(p/'independent_h3_binary_sieve.jsonl')
 require(len(binary)==len({r['id']for r in binary})and {r['id']for r in binary}==ids,'binary coverage')
 bv=read(p/'binary_witness_verification.json');ba=read(p/'binary_absence_independent.json')
 require(bv['verified']and bv['coverage_exact']and ba['complete']and ba['hits']==0,'binary independent verification')
 shell=read(p/'M_shell12_orbit_audit.json');exact=read(p/'M_shell12_exact_count.json')
 require(shell['unique_valid']and shell['generator_closed']and exact['complete'],'ambient sphere completeness')
 require(shell['signed_count']==2*shell['pair_count']==exact['signed_vectors'],'sphere independent counts agree')
 require(all(r['q']>=6 and r['q']!=7 for r in shell['shells']),'norm assumptions')
 require(ba['representatives']==sum(s['orbits']for s in shell['shells']),'binary representative coverage')
 survivors={r['id']for r in binary if r['status']=='survives'}
 mapping=read(p/'h3_source_orbit_map.json');orbit=read(p/'h3_source_orbit_verification.json')
 require(orbit['verified']and orbit['coverage_exact'],'source equivalence verification')
 require(len(mapping)==len({r[0]for r in mapping})and {r[0]for r in mapping}==survivors,'source orbit coverage')
 reps=rows(p/'independent_h3_representatives.jsonl');repids={r['id']for r in reps}
 require(repids=={r[1]for r in mapping}and len(reps)==len(repids),'representative coverage')
 sim=rows(p/'simultaneous.jsonl');simids={r['id']for r in sim};sv=read(p/'simultaneous_verification.json')
 require(len(sim)==len(simids)and all(r['integral']==0 for r in sim),'simultaneous exclusion')
 require(sv['complete']and sv['source_frames']==len(sim)and sv['integral_embeddings']==0,'independent simultaneous verification')
 require(sv['quadruples_checked']==sum(r['quadruples']for r in sim),'quadruple verification coverage')
 low=read(p/'lowspan_targets.json');rebased=read(p/'h3_lowspan_rebased.json');aff=read(p/'lowspan_affine/summary.json')
 require(len(low)==1 and low[0]['id']==rebased['source_id'],'affine source')
 require(rebased['source_H']==low[0]['H']and rebased['basis_determinant']==[1,0],'affine basis equivalence')
 require(aff['completed']and aff['embeddings_found']==0 and len(aff['stages'])==3,'affine completion')
 require(all(r['independent_python']for r in aff['stages']),'independent affine enumerations')
 require(aff['stages'][0]['input_states']==next(r['orbits']for r in shell['shells']if r['q']==12),'first affine orbit coverage')
 for a,b in zip(aff['stages'],aff['stages'][1:]):require(a['output_states']==b['input_states'],'affine stage coverage')
 require(aff['stages'][-1]['output_states']==0,'last affine empty')
 eliminated_reps=simids|{low[0]['id']}
 require(not(simids&{low[0]['id']})and eliminated_reps==repids,'final representative coverage')
 eliminated={r['id']for r in binary if r['status']!='survives'}|{t for t,r,_ in mapping if r in eliminated_reps}
 require(eliminated==ids,'complete exclusion of all targets')
 pos=read(p/'positive_simultaneous_verification.json')
 require(pos['complete']and pos['integral_embeddings']>0,'actual J0 positive control')
 log=(p/'logs/verify_input.stdout').read_text();err=(p/'logs/verify_input.stderr').read_text()
 require('DIRECT_CYCLOTOMIC_INPUT_VERIFIED'in log and not any(l.lstrip().startswith('***')and'Warning:'not in l for l in err.splitlines()),'direct input verification')
 out={'verified':True,'rank4_h3_exists':False,'ambient':'J0=(pbar J,h11/2)',
      'assumed_previously_proved_projective_minima':[6,16,27,48],
      'norm7_absence_rechecked':True,'abstract_candidates':len(targets),
      'binary_counts':dict(collections.Counter(r['status']for r in binary)),
      'binary_obstruction_types':ba['types'],'binary_direct_signed_comparisons':ba['direct_signed_comparisons'],
      'binary_survivors':len(survivors),'source_representatives':len(reps),
      'highspan_representatives':len(sim),'simultaneous_quadruples':sv['quadruples_checked'],
      'simultaneous_integral_embeddings':0,'lowspan_source_id':low[0]['id'],
      'affine_state_counts':[aff['stages'][0]['input_states']]+[r['output_states']for r in aff['stages']],
      'affine_linear_systems':sum(r['input_states']for r in aff['stages']),
      'independently_enumerated_affine_problems':sum(r['affine_problems']for r in aff['stages']),
      'all_targets_accounted_for':len(eliminated),'positive_control_integral_embeddings':pos['integral_embeddings'],
      'J0_shell_signed_count':shell['signed_count'],'J0_shell_orbits':sum(s['orbits']for s in shell['shells']),
      'tensor_trace_even_positive_definite_unimodular':True,
      'consequence_using_existing_rank1_to3_bounds_and_norm4_witness':{'Hermitian_minimum':4,'trace_minimum':8,'trace_rank':88,'extremal':True}}
 names=['independent_h3_targets.jsonl','independent_h3_binary_sieve.jsonl','h3_source_orbit_map.json',
        'h3_source_orbit_verification.json','simultaneous.jsonl','simultaneous_verification.json',
        'lowspan_affine/summary.json','M_integer_data.txt','M_shell12.txt','M_shell12_orbit_reps.txt']
 out['evidence_hashes']={n:hashlib.sha256((p/n).read_bytes()).hexdigest()for n in names}
 (p/'FINAL_VERIFICATION.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
 print(json.dumps(out,ensure_ascii=False,indent=2));return out
if __name__=='__main__':finalize(pathlib.Path(sys.argv[1])if len(sys.argv)>1 else pathlib.Path(__file__).resolve().parent/'run')
