#!/usr/bin/env python3
"""Convert fresh exact PARI exports to portable JSON and integer streams."""
import json,pathlib,sys
root=pathlib.Path(sys.argv[1])
def matrix(name):
 v=list(map(int,(root/(name+'.txt')).read_text().split()));n,m=v[:2];assert len(v)==2+n*m
 return [v[2+i*m:2+(i+1)*m]for i in range(n)]
data={'G':matrix('G88'),'M5':matrix('M5'),'M23':matrix('M23')}
(root/'lattice.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
with(root/'integer_input.txt').open('w')as f:
 f.write('88\n')
 for M in data.values():f.writelines(' '.join(map(str,r))+'\n'for r in M)
with(root/'lattice.gp').open('w')as f:
 for k,M in data.items():f.write(k+'=['+';'.join(','.join(map(str,r))for r in M)+'];\n')
