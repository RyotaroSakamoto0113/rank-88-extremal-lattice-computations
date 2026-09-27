default(nbthreads,1);
default(parisizemax,1000000000);
need(b,s)={if(!b,error(s))};
main()={my(start=getwalltime(),z,s,w,delta,bs,C,H,B,Hb,T,iso,A,autos,V,S,k,U,found=0,G7);
 z=Mod(x,polcyclo(7));s=sum(i=1,6,kronecker(i,7)*z^i);w=(1+s)/2;delta=2-z-z^-1;bs=[1,z,z^2];C=matrix(6,6,i,j,polcoef(lift(bs[(j+1)\2]*if(j%2,1,w)),i-1));need(abs(matdet(C))==1,"O_F basis of O_zeta7");
 H=matrix(3,3,i,j,sum(k=1,6,if(kronecker(k,7)==1,subst(lift(bs[i]*subst(lift(bs[j]),x,z^-1)/delta),x,z^k),0)));
 B=matrix(6,6,i,j,my(v=H[(i+1)\2,(j+1)\2]*if(i%2,1,w)*if(j%2,1,1-w));(7*polcoef(lift(v),0)-subst(lift(v),x,1))/3);
 Hb=[2,w,-1;1-w,2,w;-1,1-w,2];T=matrix(6,6,i,j,my(v=Hb[(i+1)\2,(j+1)\2]*if(i%2,1,w)*if(j%2,1,1-w));(7*polcoef(lift(v),0)-subst(lift(v),x,1))/3);
 need(denominator(B)==1 && matdet(B)==7^3 && matdet(Hb)==1,"Barnes determinants");iso=qfisom(T,B);need(type(iso)=="t_MAT" && iso~*B*iso==T,"Barnes trace isometry");
 A=matrix(6,6,i,j,if(i%2==1 && j==i+1,-2,if(i%2==0 && j==i-1,1,if(i==j && i%2==0,1,0))));
 autos=qfauto(T);V=List([matid(6)]);S=Map();mapput(S,matid(6),1);k=1;
 while(k<=#V && !found,U=iso*V[k];if(U*A==A*U,found=1;break);for(j=1,#autos[2],G7=autos[2][j]*V[k];if(!mapisdefined(S,G7),mapput(S,G7,1);listput(V,G7)));k++);
 need(found,"Barnes F-linear isometry");write(Str(OUT,"/Barnes_isometry.gp"),Str("CYCLOTOMIC_H=",H,";BARNES_H=",Hb,";CYCLOTOMIC_TRACE_G=",B,";BARNES_TRACE_G=",T,";OMEGA_ACTION=",A,";ISOMETRY=",U,";"));print("BARNES_HERMITIAN_ISOMETRY_VERIFIED");

print("BARNES_COMPLETE milliseconds=",getwalltime()-start);
};
main();
quit;
