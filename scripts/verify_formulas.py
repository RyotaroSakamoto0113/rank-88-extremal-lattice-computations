#!/usr/bin/env python3
"""Independent exact check of theta coefficients, harmonic moments and bounds."""
from fractions import Fraction as F
import json

def conv(a,b):
    return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(5)]
def power(a,n):
    r=[1,0,0,0,0]
    for _ in range(n):r=conv(r,a)
    return r
def harmonic(terms):
    out={}
    for (j,k),v in terms.items():
        if j:out[j-1,k]=out.get((j-1,k),0)+v*2*j*(88+2*j+2*k-2)
        if k>=2:out[j,k-2]=out.get((j,k-2),0)+v*8*k*(k-1)
    assert all(v==0 for v in out.values()),out
def verify():
    e=[1]+[240*sum(d**3 for d in range(1,n+1)if n%d==0)for n in range(1,5)]
    delta=[0,1,-24,252,-1472]
    basis=[conv(power(delta,j),power(e,11-3*j))for j in range(4)]
    coeff=[1,-2640,1813680,-244992000]
    theta=[sum(coeff[j]*basis[j][i]for j in range(4))for i in range(5)]
    assert theta==[1,0,0,0,168498000]
    N=theta[4];s2=F(8*N,11);s4=F(3*8**4*N,88*90)
    assert s2==122544000 and s4==261427200 and s4%9==6
    harmonic({(0,2):F(1),(1,0):F(-1,11)})
    harmonic({(0,4):F(1),(1,2):F(-12,23),(2,0):F(8,345)})
    harmonic({(0,6):F(1),(1,4):F(-5,4),(2,2):F(15,47),(3,0):F(-10,1081)})
    assert (F(960,47)*s2-F(5120,1081)*N)/10==170496000
    # Each row is a polynomial in a0,a1,a2,a3,a4 plus a constant.
    rows=[[1,2,2,2,2,2-N],[0,2,8,18,32,128-s2],
          [0,2,32,162,512,8192],[0,2,128,1458,8192,524288]]
    # S6-10*S4+1704960000=0. Eliminate a1,a2,a3.
    rel=[rows[3][i]-10*rows[2][i]for i in range(6)];rel[-1]+=1704960000
    goal=[9,0,0,0,-630,-441661950]
    goal=[goal[i]-rows[2][i] for i in range(6)]
    # 9 R0 - 49/4 R2 - (1/4)(S6 - 10 S4 + c).
    assert [9*rows[0][i]-F(49,4)*rows[1][i]-F(1,4)*rel[i]for i in range(6)]==goal
    assert 630%9==0 and 441661950%9==0
    # Tensor bounds exceed 3 for ranks 1..3, and exceed 2 for rank 4.
    products=[F(6),F(80,23),F(27,23),F(48,529)]
    assert all((s**s)*products[s-1]>3**s for s in [1,2,3])
    assert 4**4*products[3]>2**4
    assert 23**2*47<13**4  # real Hermite bound for the rank-4 first vector
    assert 46*26<11**3 and F(23*15,5)<9**2
    return {'verified':True,'theta':theta,'shell_size':N,'S2':int(s2),
            'hypothetical_S4':int(s4),'S4_residue_if_design':6,
            'harmonic_degrees_verified':[2,4,6],'tensor_bounds_verified':True}
if __name__=='__main__':print(json.dumps(verify(),indent=2))
