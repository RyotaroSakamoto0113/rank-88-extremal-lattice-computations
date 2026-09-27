#!/usr/bin/env python3
"""Bind every completed independent certificate to its actual anchor inputs."""
import argparse,csv,datetime,hashlib,json,resource,time
from pathlib import Path

NAMES=['raw_rank2','nonsaturated_rank2','ordinary_rank2_rejections','duplicate_rank2','kept_rank2','raw_rank3','nonsaturated_rank3','ordinary_rank3_rejections','second_order_rank3_rejections','duplicate_rank3','kept_rank3','raw_rank4','exact_forms']
def need(ok,s):
    if not ok:raise ValueError(s)

def main():
    p=argparse.ArgumentParser();p.add_argument('audit',type=Path);a=p.parse_args();audit=a.audit.resolve();production=audit.parent
    start=time.perf_counter();cpu=resource.getrusage(resource.RUSAGE_SELF);stamp=datetime.datetime.now(datetime.timezone.utc).isoformat();summary=json.loads((audit/'SUMMARY.json').read_text());prod=json.loads((production/'D4_SUMMARY.json').read_text())
    need(summary['status']=='verified' and summary['all78_verified'] is True and summary['verified_anchor_ids']==list(range(78)),'complete78 aggregate');need(prod['status']=='complete','production complete')
    digest_cache={}
    def sha(p):
        p=Path(p)
        if p not in digest_cache:
            h=hashlib.sha256()
            with p.open('rb') as f:
                for b in iter(lambda:f.read(1048576),b''):h.update(b)
            digest_cache[p]=h.hexdigest()
        return digest_cache[p]
    total=[0]*13;events=[];rows=[];coverage_total=0;outeruser=outersys=0
    for anchor in range(78):
        folder=audit/'anchors'/f'{anchor:03d}';actual=json.loads((folder/'verification.json').read_text());record=summary['anchors'][str(anchor)]
        need(actual['status']=='verified' and actual['anchor']==anchor,'actual verification anchor identity');need(Path(actual['shard']).resolve()==production/'anchors'/f'{anchor:03d}','actual input directory identity');need(actual['counts']==record['counts'],'per-anchor and aggregate counts agree')
        need(actual['no_rank4_candidates'] and actual['counts'][11]==0,'no rank4 candidate');need(actual['all_raw_vectors_covered_once'] and actual['all_pruning_witnesses_valid'] and actual['all_quotient_forms_reconstructed'] and actual['independent_complete_counts_match'],'all independent verification obligations')
        for path,wanted in actual['input_hashes'].items():need(sha(path)==wanted,'saved input hash unchanged: '+path)
        seen=set();stages={2:0,3:0,4:0}
        with (folder/'coverage.csv').open() as f:
            for r in csv.DictReader(f):
                key=(int(r['quotient_id']),int(r['vector_index']));need(key not in seen,'duplicate coverage row');seen.add(key);stages[int(r['rank'])]+=1
        need([stages[2],stages[3],stages[4]]==[actual['counts'][0],actual['counts'][5],actual['counts'][11]],'coverage row counts');coverage_total+=len(seen)
        total=[x+y for x,y in zip(total,actual['counts'])];d=record['dispatch'];gp=actual['gp_stage'];outeruser+=d['wait4_user_cpu_seconds'];outersys+=d['wait4_system_cpu_seconds'];events.extend([(datetime.datetime.fromisoformat(d['started_utc']),1),(datetime.datetime.fromisoformat(d['finished_utc']),-1)])
        rows.append({'anchor_id':anchor,'started_utc':d['started_utc'],'finished_utc':d['finished_utc'],'dispatch_wall_seconds':d['wall_seconds'],'verification_wall_seconds':actual['total_wall_seconds'],'gp_wall_seconds':gp['wall_seconds'],'gp_user_cpu_seconds':gp['user_cpu_seconds'],'gp_system_cpu_seconds':gp['system_cpu_seconds'],'python_user_cpu_seconds':actual['python_user_cpu_seconds'],'python_system_cpu_seconds':actual['python_system_cpu_seconds'],'wait4_user_cpu_seconds_including_startup':d['wait4_user_cpu_seconds'],'wait4_system_cpu_seconds_including_startup':d['wait4_system_cpu_seconds'],'rank2_raw':actual['counts'][0],'rank3_raw':actual['counts'][5],'rank4_raw':actual['counts'][11],'exact_forms':actual['counts'][12]})
    need(total==summary['counts']==[prod['counts'][name] for name in NAMES],'aggregate counts bind to production')
    active=peak=0
    for t,d in sorted(events):active+=d;peak=max(peak,active);need(active>=0,'valid concurrency intervals')
    need(active==0 and peak<=4,'bounded4-way concurrency')
    with (audit/'per_anchor_timing.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    endcpu=resource.getrusage(resource.RUSAGE_SELF);result={'status':'verified','actual_anchor_ids':list(range(78)),'all_actual_anchor_ids_match_directories':True,'all_input_hashes_unchanged':True,'unique_input_files_hashed':len(digest_cache),'all_coverage_rows_unique':True,'coverage_rows':coverage_total,'aggregate_counts_match_production':True,'counts':dict(zip(NAMES,total)),'peak_concurrent_verifiers':peak,'collector_elapsed_wall_seconds':summary['elapsed_wall_including_waits_seconds'],'sum_gp_and_python_user_cpu_seconds_excluding_interpreter_startup':summary['sum_user_cpu_seconds'],'sum_gp_and_python_system_cpu_seconds_excluding_interpreter_startup':summary['sum_system_cpu_seconds'],'sum_wait4_user_cpu_seconds_including_startup':outeruser,'sum_wait4_system_cpu_seconds_including_startup':outersys,'timing':{'started_utc':stamp,'wall_seconds':time.perf_counter()-start,'user_cpu_seconds':endcpu.ru_utime-cpu.ru_utime,'system_cpu_seconds':endcpu.ru_stime-cpu.ru_stime},'audit_source_sha256':sha(Path(__file__))}
    (audit/'final_binding_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
