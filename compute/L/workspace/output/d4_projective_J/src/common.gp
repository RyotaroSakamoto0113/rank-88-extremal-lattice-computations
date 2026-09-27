\\ Exact quotient algebra shared by the new independent minima audit.
need(b,s)={if(!b,error(s))};
read(Str(OUT,"/ambient.gp"));read(Str(OUT,"/generators.gp"));
n=22;Id=matid(n);NN=[1,2,2];
need(G==G~ && denominator(G)==1 && matdet(G)==23^11,"ambient trace Gram");
need(W^2-W+6*Id==0 && W~*G==G*(Id-W),"omega adjoint");
need(denominator((G-2*W~*G)/23)==1 && denominator((11*G+W~*G)/23)==1,"O-integrality");
for(i=1,#AUT[2],A=AUT[2][i];need(denominator(A)==1 && abs(matdet(A))==1 && A~*G*A==G && A*W==W*A,"generator is O-isometry"));
twist(R,c)={if(c==1,return(R));mathnf(concat(2*R,(W-if(c==2,Id,0))*R))/2};
CC=vector(3,c,twist(Id,c));GG=vector(3);CI=vector(3);
for(c=1,3,L=NN[c]*CC[c]~*G*CC[c];U=qflllgram(L);need(abs(matdet(U))==1,"cusp LLL unimodular");CC[c]=CC[c]*U;CI[c]=CC[c]^-1;GG[c]=NN[c]*CC[c]~*G*CC[c];need(denominator(GG[c])==1 && GG[c]==GG[c]~,"cusp Gram integral"));
isprimitive(K)={my(r=matsize(K)[2],sn);need(denominator(K)==1,"noninteger module");sn=matsnf(K);sum(i=1,#sn,sn[i]!=0)==r && vecmax(abs(sn))==1};
delta(K)={my(r=matsize(K)[2]/2,D=matdet(K~*G*K)/23^r,d);need(denominator(D)==1 && D>0,"determinant square not positive integer");d=sqrtint(D);need(d*d==D,"determinant square mismatch");d};
linebasis(y,c)={if(c==1,return(Mat([y,W*y])));if(c==2,Mat([2*y,W*y]),Mat([2*y,(W-Id)*y]))};
schur(R)={G-G*R*(R~*G*R)^-1*R~*G};
quotient(R,rho,c)={
 my(K=CI[c]*twist(R,c),k=matsize(R)[2],sn,U,T,Q,Qr,L,den);
 need(denominator(K)==1 && isprimitive(K),"twisted kernel must be primitive");
 sn=matsnf(K,1);U=sn[1]^-1;T=U[,1..n-k];
 need(abs(matdet(concat(K,T)))==1,"quotient complement not unimodular");
 Q=rho*NN[c]*CC[c]~*schur(R)*CC[c];need(Q*K==0,"wrong quotient kernel");
 Qr=T~*Q*T;L=qflllgram(Qr);need(abs(matdet(L))==1,"quotient LLL unimodular");T=T*L;Qr=L~*Qr*L;
 den=denominator(Qr);Qr*=den;need(Qr==Qr~ && denominator(Qr)==1,"quotient not integral");
 [Qr,T,den]
};
canon(v)={for(i=1,#v,if(v[i],if(v[i]<0,return(-v),return(v))));v};
validate(Q,b,sv)={my(V=sv[3],S);need(sv[1]==2*matsize(V)[2],"candidate count");S=Set(vector(matsize(V)[2],i,canon(V[,i])));need(#S==matsize(V)[2],"duplicate sign pairs");for(i=1,#S,need(S[i]~*Q*S[i]>0 && S[i]~*Q*S[i]<=b,"invalid vector"));1};
batchno=0;fdforms=fileopen(Str(OUT,"/exact_inputs.txt"),"w");fdmanifest=fileopen(Str(OUT,"/exact_manifest.csv"),"w");fdvectors=fileopen(Str(OUT,"/enumerated_vectors.gp"),"w");
emit(Q,b,sv,meta)={my(nq=matsize(Q)[1]);validate(Q,b,sv);batchno++;filewrite(fdforms,Str(nq," ",b));for(i=1,nq,filewrite(fdforms,strjoin(vector(nq,j,Str(Q[i,j]))," ")));filewrite(fdmanifest,Str(batchno,",",nq,",",b,",",sv[1],",",meta));filewrite(fdvectors,Str("V",batchno,"=",sv[3],";"));batchno};
closecert()={fileclose(fdforms);fileclose(fdmanifest);fileclose(fdvectors)};
