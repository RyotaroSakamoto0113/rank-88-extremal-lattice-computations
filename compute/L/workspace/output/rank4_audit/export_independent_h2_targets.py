import json,pathlib
p=pathlib.Path(__file__).parent
source=(p/'independent_h2_abstract.py').read_text().split('start=time.monotonic()')[0]
scope={'__file__':str(p/'independent_h2_abstract.py')}
exec(compile(source,str(p/'independent_h2_abstract.py'),'exec'),scope)
mul,neg,transform=map(scope.get,('mul','neg','transform'))
d=scope['D'];zero=(0,0)
E=[[mul(neg(d),x)for x in r]for r in (((1,0),zero,(4,0),(1,-1)),(zero,(1,0),(0,-1),(4,0)))]+[[zero,zero,(23,0),zero],[zero,zero,zero,(23,0)]]
out=[]
for i,c in enumerate(json.loads((p/'independent_h2_abstract.json').read_text())['candidates']):
 HX=transform(E,c['HY'])
 assert all(a%529==b%529==0 for r in HX for a,b in r)
 HX=[[[a//529,b//529]for a,b in r]for r in HX]
 out.append(dict(id=i,H=HX,det=c['detX']))
(p/'independent_h2_targets.jsonl').write_text('\n'.join(json.dumps(c)for c in out)+'\n')
print('Exported',len(out),'integral Hermitian targets.')
