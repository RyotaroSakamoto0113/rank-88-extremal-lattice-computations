"""Reproduce every J0 exclusion in a fresh directory, from source and inputs."""
import argparse,datetime,json,pathlib,shutil,sys,tempfile,time
import run_initial as pipeline
from finalize import finalize
HERE=pathlib.Path(__file__).resolve().parent
SOURCES=['independent_chart_audit.gp','export_charts_json.gp','ambient_qfauto.gp',
 'check_M_shell_orbits.cpp','check_binary_absence.cpp','independent_embedding_sieve.cpp',
 'independent_simultaneous.cpp','verify_simultaneous.py','h3_source_orbits.gp',
 'exact_affine_enum.cpp','affine_embedding_audit.py','affine_metadata_verify.py',
 'affine_stage_generator.gp','test_affine_enum.py','abstract_h3.cpp',
 'verify_binary_witnesses.py','verify_source_orbits.py','prepare_lowspan_target.py','lowspan_seed.json']
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run-dir',type=pathlib.Path);args=ap.parse_args()
 if args.run_dir:root=args.run_dir.resolve();root.mkdir(parents=True,exist_ok=False)
 else:
  base=HERE/'reproductions';base.mkdir(exist_ok=True);root=pathlib.Path(tempfile.mkdtemp(prefix=datetime.datetime.now().strftime('%Y%m%d-%H%M%S-'),dir=base))
 out=root/'output/rank4_audit';out.mkdir(parents=True);(out/'logs').mkdir()
 (root/'input').mkdir();(root/'output/d3/benchmark').mkdir(parents=True)
 for n in SOURCES:shutil.copy2(HERE/'run/output/rank4_audit'/n,out/n)
 for n in ['J_ambient.gp','J0_exact_input.gp','norm4_witness.json']:shutil.copy2(HERE/'run/input'/n,root/'input'/n)
 shutil.copy2(HERE/'run/output/d3/benchmark/exact_enum.cpp',root/'output/d3/benchmark/exact_enum.cpp')
 pipeline.ROOT=root;pipeline.OUT=out;pipeline.timings=[]
 started=time.perf_counter();print('FRESH_RUN',root,flush=True)
 pipeline.main()
 run=pipeline.run
 run('verify_input',['gp','-fq',HERE/'verify_input.gp'])
 run('prepare_lowspan',[sys.executable,out/'prepare_lowspan_target.py'])
 run('affine_enumeration',[sys.executable,out/'affine_embedding_audit.py',out/'h3_lowspan_rebased.json','lowspan_affine'])
 run('affine_regression',[sys.executable,out/'test_affine_enum.py'])
 w=json.loads((root/'input/norm4_witness.json').read_text());h=[[[a,b]for a,b in zip(ra,rb)]for ra,rb in zip(w['HX_a'],w['HX_b'])]
 control={'id':-1,'det':96,'H':h}
 (out/'positive_control.jsonl').write_text(json.dumps(control)+'\n');pipeline.flat([control],out/'positive_control.flat')
 run('positive_binary',[out/'independent_embedding_sieve',out/'positive_control.flat'],out/'positive_binary.jsonl')
 run('positive_simultaneous',[out/'independent_simultaneous',out/'positive_control.flat'],out/'positive_simultaneous.jsonl')
 run('positive_verify',[sys.executable,out/'verify_simultaneous.py','--targets',out/'positive_control.jsonl','--results',out/'positive_simultaneous.jsonl','--quadruples',str(out/'positive_control.flat')+'.baseline_quads','--output',out/'positive_simultaneous_verification.json'])
 result=finalize(root);result['fresh_run_wall_seconds']=time.perf_counter()-started
 (out/'FULL_RUN_SUMMARY.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
 print('PROOF_COMPLETE',out/'FULL_RUN_SUMMARY.json',flush=True)
if __name__=='__main__':main()
