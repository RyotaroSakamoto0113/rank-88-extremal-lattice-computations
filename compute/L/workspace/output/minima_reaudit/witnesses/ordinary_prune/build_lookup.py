#!/usr/bin/env python3
"""Pack a completely known sphere orbit map into collision-free GP integers.

This is a data conversion audit, not a fresh sphere enumeration. It checks
that the packed lookup contains exactly the provided full sphere, modulo sign.
"""
import argparse,datetime,hashlib,json,resource,time
from pathlib import Path

def need(ok,s):
    if not ok:raise ValueError(s)

def canon(v):
    first=next((x for x in v if x),0);need(first!=0,'nonzero vector')
    return tuple(v) if first>0 else tuple(-x for x in v)

def pack(v):
    need(len(v)==22 and all(-32<=x<32 for x in v),'base64 coordinate guard')
    key=0
    for x in v:key=64*key+x+32
    return key

def unpack(k):
    v=[]
    for i in range(22):v.append(k%64-32);k//=64
    need(k==0,'packed key width');return tuple(reversed(v))

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    root=Path(__file__).resolve().parent;project=root.parents[3]
    a=argparse.ArgumentParser();a.add_argument('--lookup',type=Path,default=root.parent/'orbits/point_orbit_lookup.txt');a.add_argument('--ambient',type=Path,default=root.parent/'M_integer_data.txt');a.add_argument('--sphere',type=Path,default=project/'output/rank4_d27_d48/runs/20260912-175210-2mzyf0pb/output/rank4_audit/M_shell12.txt');a.add_argument('--output',type=Path,default=root);p=a.parse_args();p.output.mkdir(parents=True,exist_ok=True)
    started=datetime.datetime.now(datetime.timezone.utc).isoformat();t=time.perf_counter();c=resource.getrusage(resource.RUSAGE_SELF)
    items={};raw=[];coords=[];norms={};orbits={}
    with p.lookup.open() as f:
        nv,n=map(int,next(f).split());need(n==22,'lookup dimension')
        for line in f:
            r=list(map(int,line.split()));need(len(r)==24,'lookup row width');q,orb=r[:2];v=tuple(r[2:]);need(v==canon(v),'lookup canonical sign');need(1<=q<=12 and 0<=orb<128,'packed value guard');k=pack(v);need(unpack(k)==v,'exact pack/unpack');need(k not in items,'injective packed keys');items[k]=128*q+orb;raw.append(v);coords.extend(v);norms[q]=norms.get(q,0)+1;orbits.setdefault(orb,q);need(orbits[orb]==q,'orbit has unique norm')
    need(len(items)==nv,'complete lookup rows');known=set()
    with p.sphere.open() as f:
        ns,n=map(int,next(f).split());need(n==22 and ns==nv,'same sphere count')
        for line in f:
            v=canon(tuple(map(int,line.split())));k=pack(v);need(k not in known,'sphere unique modulo sign');known.add(k)
    need(len(known)==ns and known==set(items),'lookup key set equals entire known sphere')
    data=list(map(int,p.ambient.read_text().split()));need(data[0]==22,'ambient rank');G=[data[1+22*i:1+22*(i+1)] for i in range(22)]
    entries=sorted(items.items());datafile=p.output/'ordinary_lookup.gp'
    with datafile.open('w') as f:
        f.write('ORDINARY_G=['+';'.join(','.join(map(str,row)) for row in G)+'];\n')
        f.write('ORDINARY_KEYS=['+','.join(str(k) for k,v in entries)+'];\n')
        f.write('ORDINARY_VALUES=['+','.join(str(v) for k,v in entries)+'];\n')
        f.write('ORDINARY_MAX_NORM=12;\nORDINARY_PAIR_COUNT='+str(nv)+';\nORDINARY_ORBIT_COUNT='+str(len(orbits))+';\n')
    d=resource.getrusage(resource.RUSAGE_SELF);meta={'status':'verified_conversion','not_a_new_sphere_enumeration':True,'base':64,'coordinate_offset':32,'coordinate_guard':[-32,31],'observed_coordinate_range':[min(coords),max(coords)],'count':nv,'orbit_count':len(orbits),'encoded_key_bits':max(k.bit_length() for k in items),'injective':True,'pack_unpack_identity':True,'same_keys_as_entire_known_sphere':True,'shell_counts':norms,'inputs':{str(x):sha(x) for x in [p.lookup,p.ambient,p.sphere]},'output':{'path':str(datafile),'sha256':sha(datafile),'bytes':datafile.stat().st_size},'timing':{'started_utc':started,'wall_seconds':time.perf_counter()-t,'user_cpu_seconds':d.ru_utime-c.ru_utime,'system_cpu_seconds':d.ru_stime-c.ru_stime}}
    (p.output/'lookup_conversion.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))

if __name__=='__main__':main()
