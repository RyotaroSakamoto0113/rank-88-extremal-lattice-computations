"""Join independently verified stages, rejecting missing or duplicate targets."""
from pathlib import Path
import collections, hashlib, json
P=Path(__file__).resolve().parent
def read(name):return json.loads((P/name).read_text())
def lines(name):return [json.loads(s)for s in(P/name).read_text().splitlines()]
def unique(rows):
 d={r['id']:r for r in rows};assert len(d)==len(rows);return d
def main():
 targets=unique(lines('independent_h3_targets.jsonl'));cheap=unique(lines('independent_h3_binary_sieve.jsonl'))
 assert len(targets)==13272 and set(targets)==set(cheap)
 c=read('binary_witness_verification.json');assert c['verified'] and c['candidates']==13272
 a=read('binary_absence_independent.json');assert a['complete'] and a['hits']==0 and a['types']==c['distinct_binary_types']==126
 high={i for i,r in cheap.items()if r['status']=='survives' and r['short_rank']==4}
 low={i for i,r in cheap.items()if r['status']=='survives' and r['short_rank']<4}
 assert len(high)==1365 and len(low)==12
 sim=unique(lines('independent_h3_simultaneous_baseline.jsonl'))
 assert all(r['integral']==0 for r in sim.values())
 verification_name='independent_h3_simultaneous_verification.json' if (P/'independent_h3_simultaneous_verification.json').exists() else 'simultaneous_independent_verification.json'
 sv=read(verification_name)
 assert sv['complete'] and sv['source_frames']==len(sim) and sv['quadruples_checked']==sum(r['quadruples']for r in sim.values())
 assert sv['integral_embeddings']==0 and sv['all_gram_conditions_match'] and sv['field_inverse_matches_restriction_inverse'] and sv['positive_integral_reconstruction_controls']==len(sim)
 orbit=read('h3_source_orbit_map.json');om={i:(rep,k)for i,rep,k in orbit}
 assert len(om)==len(orbit)==1377 and set(om)==high|low
 ov=read('h3_source_orbit_verification.json');assert ov['verified'] and ov['targets']==1377 and ov['representatives']==363
 for name,expected in ov['hashes'].items():assert hashlib.sha256((P/name).read_bytes()).hexdigest()==expected
 high_reps={om[i][0] for i in high}
 assert len(high_reps)==362 and set(sim) in (high,high_reps)
 assert {om[i][0]for i in low}=={11622}
 L=read('lowspan_affine/summary.json');assert L['completed'] and L['embeddings_found']==0
 assert [r['input_states']for r in L['stages']]==[45,49,9]
 assert [r['output_states']for r in L['stages']]==[49,9,0]
 assert sum(r['affine_problems']for r in L['stages'])==101 and all(r['independent_python']for r in L['stages'])
 H2=unique(lines('independent_h2_targets.jsonl'));hw=read('h2_candidate_witnesses_summary.json')
 assert len(H2)==23 and hw['candidates']==23 and hw['all_excluded']
 hm=read('h2_M_result.json');hn=read('h2_N_result.json')
 assert hm['complete'] and hm['pairs']==10626 and hm['obstruction_representations']==[0,0,0,0]
 assert hn['complete'] and hn['obstruction_representations']==[0,0,6,6]
 if (P/'d3_dependency_summary.json').exists():
  d3=read('d3_dependency_summary.json');assert d3['complete'] and d3['forms']==30423 and d3['rank18_vectors']==0
  dependency_forms=30423
 else:
  dep=read('FULL_RUN_SUMMARY.json')['dependency']
  assert dep['input_sha256']=='eed37d27b0950a5c97bd03f59767e116545daf2cec4baaa234c8c3a84e92dc8e'
  dependency_forms=0
  if (P/'d3_dependency_recheck.jsonl').exists():
   dd=lines('d3_dependency_recheck.jsonl');assert len(dd)==30423 and all(r['complete'] for r in dd)
   assert all(r['signed_vectors']==0 for r in dd if r['dimension']==18)
   dependency_forms=30423
 ambient_log=P/'ambient_input_verification.log'
 if not ambient_log.exists():ambient_log=P/'driver_logs/verify_ambient_input.stdout'
 assert 'AMBIENT_INPUT_VERIFIED' in ambient_log.read_text()
 if (P/'d3_shell_reaudit.log').exists():assert 'ALL_SAVED_SHELL_ORBITS_VERIFIED' in(P/'d3_shell_reaudit.log').read_text()
 if (P/'d3_quotients18_reaudit.log').exists():assert 'min_det_R=16' in(P/'d3_quotients18_reaudit.log').read_text()
 sphere=read('M_shell12_orbit_audit.json');assert sphere['unique_valid'] and sphere['generator_closed'] and sphere['signed_count']==567226
 exact=read('M_shell12_exact_count.json');assert exact['complete'] and exact['signed_vectors']==567226
 summary={'claim':'For the specified M and I, no F-rank-4 tensor has h(z,z)<=3.','verified':True,
  'h1':'excluded by integrality and 23h>=24','h2_candidates':23,'h2_remaining':0,'h3_candidates':len(targets),
  'h3_cheap_exclusions':dict(collections.Counter(r['status']for r in cheap.values())),
  'h3_source_representatives':363,'h3_simultaneous_targets_checked':len(sim),'h3_simultaneous_original_targets_covered':len(high),'h3_lowspan_targets':len(low),
  'h3_lowspan_representative':11622,'h3_remaining':0,'d2_d3_dependency_forms_rechecked':dependency_forms,
  'scope_note':'The independent reconstruction certifies the mathematical conclusion; it does not claim a rerun of the unavailable old ZIP, or a new isometry to a separately stored rank-88 Gram matrix.'}
 files=['independent_h2_targets.jsonl','independent_h3_targets.jsonl','independent_h3_binary_sieve.jsonl','independent_h3_simultaneous_baseline.jsonl',verification_name,'h3_source_orbit_map.json','lowspan_affine/summary.json','binary_witness_verification.json','binary_absence_independent.json','M_integer_data.txt','M_shell12_exact_count.json']
 summary['sha256']={n:hashlib.sha256((P/n).read_bytes()).hexdigest()for n in files}
 (P/'FINAL_AUDIT_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
 print(json.dumps({k:v for k,v in summary.items()if k!='sha256'},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
