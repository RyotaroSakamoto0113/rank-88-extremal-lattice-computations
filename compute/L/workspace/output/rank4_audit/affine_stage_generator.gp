/* Exact affine problems for extending prescribed Hermitian Gram embeddings. */
read("output/d3/benchmark/ambient.gp");
w=Mod(x,x^2-x+6);
assert(b,s)={if(!b,print("FAIL: ",s);quit(1))};
jq(q)=Str("\"",q,"\"");
jvec(v)=Str("[",strjoin(vector(#v,i,jq(v[i])),","),"]");
jmat(M)=Str("[",strjoin(vector(matsize(M)[1],i,jvec(Vec(M[i,]))),","),"]");
ab(h)=[polcoef(lift(h),0),polcoef(lift(h),1)];
hval(u,v)={my(t=u~*G*v,s=u~*G*W*v);(12*t-s+(2*s-t)*w)/23};
affine_stage(States,Target,prefix)={
  my(n=22,fd=fileopen(Str(prefix,"_input.txt"),"w"),fm=fileopen(Str(prefix,"_metadata.jsonl"),"w"),V,k,m,A,rhs,hn,T,H,K,B,Q,U,y,x0,c,r,shift,base,entries,tq,affid=0,consistent=0);
  for(idx=1,#States,
    V=States[idx];k=matsize(V)[2];m=n-2*k;
    assert(k>=1 && k<=3 && matsize(V)[1]==n,"state dimensions");
    for(i=1,k,for(j=1,k,assert(hval(V[,i],V[,j])==Target[i,j],"state Gram")));
    A=matrix(2*k,n,i,j,if(i<=k,(V[,i]~*G)[j],(V[,i-k]~*G*W)[j]));
    rhs=vector(2*k,i,my(z=ab(Target[if(i<=k,i,i-k),k+1]));if(i<=k,2*z[1]+z[2],z[1]+12*z[2]))~;
    hn=mathnf(A,1);H=hn[1];T=hn[2];
    assert(matsize(H)==[2*k,2*k] && abs(matdet(T))==1,"full HNF kernel");
    assert(A*T==concat(matrix(2*k,m),H),"HNF identity");
    y=H^-1*rhs;
    base=Str("{\"state_id\":",idx,",\"k\":",k,",\"A\":",jmat(A),",\"rhs\":",jvec(rhs),",\"T\":",jmat(T),",\"H\":",jmat(H));
    if(denominator(y)!=1,
      filewrite(fm,Str(base,",\"consistent\":false}"));next);
    consistent++;
    K=matrix(n,m,i,j,T[i,j]);
    B=matrix(n,2*k,i,j,T[i,j+m]);
    x0=B*y;
    Q=K~*G*K;U=qflllgram(Q);
    assert(abs(matdet(U))==1,"kernel LLL unimodular");
    K=K*U;Q=K~*G*K;
    assert(abs(matdet(concat(K,B)))==1 && A*K==0 && A*x0==rhs,"whole integer solution lattice");
    c=-Q^-1*K~*G*x0;
    tq=2*polcoef(lift(Target[k+1,k+1]),0);
    r=tq-x0~*G*x0+c~*Q*c;
    shift=vector(m,i,round(c[i]))~;
    x0=x0+K*shift;c=c-shift;
    assert(Q*c==-K~*G*x0 && r==tq-x0~*G*x0+c~*Q*c,"affine square identity");
    affid++;
    filewrite(fd,Str(m," ",r));filewrite(fd,strjoin(vector(m,i,Str(c[i]))," "));
    for(i=1,m,filewrite(fd,strjoin(vector(m,j,Str(Q[i,j]))," ")));
    entries=Str(",\"consistent\":true,\"affine_id\":",affid,",\"B\":",jmat(B),",\"K\":",jmat(K),",\"x0\":",jvec(x0),",\"Q\":",jmat(Q),",\"center\":",jvec(c),",\"radius\":",jq(r),",\"target_trace_norm\":",jq(tq));
    filewrite(fm,Str(base,entries,"}"));
  );
  fileclose(fd);fileclose(fm);
  print("states=",#States," consistent=",consistent," affine=",affid);
};
