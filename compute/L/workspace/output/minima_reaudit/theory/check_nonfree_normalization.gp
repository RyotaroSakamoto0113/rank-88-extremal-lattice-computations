default(parisize,100000000);
assert(b,s)={if(!b,error(s))};
w=Mod(u,u^2-u+6);d=2*w-1;
bar(z)=subst(lift(z),u,1-w);
tr(z)={my(f=lift(z));2*polcoef(f,0)+polcoef(f,1)};
H=[1,1/2+7*d/46,1/2+d/46;1/2-7*d/46,1,1/4-3*d/92;1/2-d/46,1/4+3*d/92,1/2];
assert(matdet(H)==1/46,"rawdet");
for(orient=0,1,B=[1,w,0,0,0,0;0,0,1,w,0,0;0,0,0,0,2,w-orient];G=matrix(6,6,i,j,tr((B[,i]~*H*bar(B[,j]))));W0=[0,-6;1,1];W1=[orient,-(orient^2-orient+6)/2;2,1-orient];W=matrix(6,6,i,j,if((i-1)\2!=(j-1)\2,0,if(i<=4,W0[(i-1)%2+1,(j-1)%2+1],W1[(i-1)%2+1,(j-1)%2+1])));assert(W^2-W+6*matid(6)==0,"Wpoly");assert(W~*G==G*(1-W),"adjoint");CC=[matid(6),mathnf(concat(2*matid(6),W-matid(6)))/2,mathnf(concat(2*matid(6),W))/2];NN=[1,2,2];for(c=1,3,Q=NN[c]*CC[c]~*G*CC[c];den=denominator(Q);Q*=den;ev=qfminim(Q,,,2);print("orient=",orient," cusp=",c," halfnorm=",ev[2]/(2*den)," vectors=",ev[1]," normalized determinant=",2*matdet(H))));
quit;
