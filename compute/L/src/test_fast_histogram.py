from pathlib import Path
import subprocess,json,random,itertools,math,time
O=Path(__file__).resolve().parents[1];R=O/'results/fast_regression';R.mkdir(exist_ok=True);rng=random.Random(880023);t=time.perf_counter()
forms=[];expected=[]
for case in range(32):
 n=case%4+1;A=[[rng.randrange(-2,3) for j in range(n)] for i in range(n)];G=[[sum(A[k][i]*A[k][j] for k in range(n))+(i==j) for j in range(n)] for i in range(n)];b=rng.randrange(1,21);h=[0]*(b+1)
 for v in itertools.product(range(-math.isqrt(b),math.isqrt(b)+1),repeat=n):
  q=sum(v[i]*G[i][j]*v[j] for i in range(n) for j in range(n))
  if 0<q<=b:h[q]+=1
 forms.append(str(n)+' '+str(b)+'\n'+'\n'.join(' '.join(map(str,r)) for r in G)+'\n');expected.append(h)
(R/'forms.txt').write_text(''.join(forms));p=subprocess.run([str(O/'results/fast_histogram'),'--batch',str(R/'forms.txt')],capture_output=True,text=True);assert p.returncode==0,p.stderr
rows=[json.loads(s) for s in p.stdout.splitlines()];assert len(rows)==32 and all(r['complete'] and r['histogram']==h for r,h in zip(rows,expected));(R/'results.jsonl').write_text(p.stdout)
for name,s in [('asymmetric','2 4\n2 1\n0 2\n'),('indefinite','2 4\n1 2\n2 1\n'),('oversize','4 42\n1000000000 0 0 0\n0 1000000000 0 0\n0 0 1000000000 0\n0 0 0 1000000000\n')]:
 (R/(name+'.txt')).write_text(s);p=subprocess.run([str(O/'results/fast_histogram'),str(R/(name+'.txt'))],capture_output=True,text=True);assert p.returncode!=0 and not p.stdout;(R/(name+'.stderr')).write_text(p.stderr)
small=[json.loads(s) for s in (O/'results/fast_nonidentity.jsonl').read_text().splitlines()];slow=[json.loads(s) for s in (O/'results/orbit_exact_remaining/nonidentity.stdout').read_text().splitlines()];assert len(small)==len(slow)==42 and all(a['histogram']==b['histogram'] for a,b in zip(small,slow))
r={'status':'verified','brute_force_forms':32,'GMP_rational_cross_checks':42,'invalid_or_unsupported_inputs_rejected':3,'wall_seconds':time.perf_counter()-t};(R/'verification.json').write_text(json.dumps(r,indent=2)+'\n');print(r)
