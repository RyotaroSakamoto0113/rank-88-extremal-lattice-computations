\\ Independent saved-certificate verification. No common.gp or pruning search
\\ function is loaded. Bareiss/cofactor arithmetic reconstructs every form.
vb_need(b,s)={if(!b,error(Str("D4 certificate: ",s)))};
vb_det(A)={
 my(B=A,n=matsize(A)[1],previous=1,sg=1,p,t,z);
 vb_need(matsize(A)[2]==n&&denominator(A)==1,"Bareiss integer square input");
 if(n==0,return(1));
 for(k=1,n-1,
  if(B[k,k]==0,p=0;for(i=k+1,n,if(!p&&B[i,k],p=i));if(!p,return(0));t=B[k,];B[k,]=B[p,];B[p,]=t;sg=-sg);
  for(i=k+1,n,for(j=k+1,n,z=(B[i,j]*B[k,k]-B[i,k]*B[k,j])/previous;vb_need(denominator(z)==1,"Bareiss exact division");B[i,j]=z));
  previous=B[k,k];for(i=k+1,n,B[i,k]=0);
 );sg*B[n,n]
};
vb_minor(A,r,c)={my(n=matsize(A)[1]);matrix(n-1,n-1,i,j,A[i+(i>=r),j+(j>=c)])};
vb_adj(A)={my(n=matsize(A)[1]);matrix(n,n,i,j,(-1)^(i+j)*vb_det(vb_minor(A,j,i)))};
vb_primitive(A)={my(s=matsnf(A),r=matsize(A)[2]);denominator(A)==1&&sum(i=1,#s,s[i]!=0)==r&&vecmax(abs(s))==1};
vb_equal_module(A,B)={
 my(U=matinverseimage(A,B));
 if(type(U)!="t_MAT"||matsize(U)!=[matsize(A)[2],matsize(B)[2]],return(0));
 A*U==B&&denominator(U)==1&&matsize(U)[1]==matsize(U)[2]&&abs(vb_det(U))==1
};
vb_in_module(A,B)={my(U=matinverseimage(A,B));type(U)=="t_MAT"&&matsize(U)==[matsize(A)[2],matsize(B)[2]]&&A*U==B&&denominator(U)==1};
vb_delta(R)={
 my(k=matsize(R)[2],d,D);
 vb_need(matsize(R)[1]==22&&k%2==0&&denominator(R)==1,"integral even-rank module");
 d=vb_det(R~*SAVED_G*R);vb_need(d>0,"module full rank");D=d/23^(k/2);vb_need(denominator(D)==1,"trace determinant normalization");d=sqrtint(D);vb_need(d*d==D,"Hermitian determinant square");d
};
vb_twist(R,c)={
 my(k=matsize(R)[2],A,K,H,B);
 if(c==1,return(R));
 A=matinverseimage(R,SAVED_W*R);vb_need(type(A)=="t_MAT"&&matsize(A)==[k,k]&&R*A==SAVED_W*R&&denominator(A)==1,"O stability for ideal inverse");
 if(c==3,A-=matid(k));
 \\ y=R*t/2 lies in a_c^-1 R exactly when A*t is divisible by2.
 K=matkerint(concat(A,-2*matid(k)));vb_need(matsize(K)==[2*k,k]&&vb_primitive(K),"primitive integer congruence kernel");
 vb_need(concat(A,-2*matid(k))*K==0,"congruence kernel identity");B=K[1..k,];
 vb_need(abs(vb_det(B))==2^matrank(Mod(A,2)),"full mod2 congruence lattice");
 R*B/2
};
vb_line(y,c)={if(c==1,Mat([y,SAVED_W*y]),if(c==2,Mat([2*y,SAVED_W*y]),Mat([2*y,(SAVED_W-matid(22))*y])))};
vb_canonical(R)={my(best=mathnf(R),H);for(i=1,#SAVED_STAB,H=mathnf(SAVED_STAB[i]*R);if(cmp(H,best)<0,best=H));best};
vb_ordinary(R,w,m,anchor)={
 my(v=w[3],z=w[4],q=w[1],oid=w[2],s=0,key,j,val);
 vb_need(#w==4&&R*z==v&&denominator(v)==1&&denominator(z)==1,"ordinary witness lies in candidate");
 vb_need((v~*SAVED_G*v)/2==q&&type(q)=="t_INT"&&q>0&&q<=12,"ordinary exact norm");
 for(i=1,22,if(!s&&v[i],s=sign(v[i])));vb_need(s!=0,"nonzero ordinary witness");
 for(i=1,22,vb_need(s*v[i]>=-32&&s*v[i]<32,"packed coordinate guard"));
 key=sum(i=1,22,(s*v[i]+32)*2^(6*(22-i)));j=vecsearch(ORDINARY_KEYS,key);vb_need(j>0,"ordinary membership in complete sphere");val=ORDINARY_VALUES[j];
 vb_need(val\128==q&&val%128==oid,"ordinary orbit identity");vb_need(q<m||(q==m&&oid<anchor),"strictly earlier ordinary vector");1
};
vb_second(P,R,T,dr,w)={
 my(c=w[1],value=w[2],y=w[3],Rp=w[4],canp=w[5],expected);
 vb_need(c>=1&&c<=3,"second-order ideal class");expected=concat(P,vb_line(y,c));
 vb_need(Rp==expected&&denominator(Rp)==1&&vb_in_module(T,Rp),"second-order witness Rp inside T and containing P");
 vb_need(vb_delta(Rp)==value&&value>0&&value<=dr,"second-order determinant");
 if(value<dr,return(1));
 vb_need(vb_primitive(Rp),"equal-determinant second module primitive");vb_need(canp==vb_canonical(Rp)&&cmp(canp,R)<0,"earlier canonical second module");1
};
vb_quotient(qid,meta,form,bound)={
 my(R=meta[1],rho=meta[2],c=meta[3],T=meta[4],den=meta[5],A,detA,adjA,U,expected,K,tw,k=matsize(R)[2]);
 vb_need(c>=1&&c<=3&&rho==vb_delta(R)&&vb_primitive(R),"quotient prefix determinant and saturation");
 vb_need(matsize(T)==[22,22-k]&&denominator(T)==1&&type(den)=="t_INT"&&den>0,"quotient complement dimensions");
 tw=vb_twist(R,c);K=VB_CI[c]*tw;vb_need(denominator(K)==1&&vb_primitive(K),"exact twisted integer kernel");
 vb_need(abs(vb_det(concat(K,T)))==1,"kernel plus complement unimodular");
 A=R~*SAVED_G*R;detA=vb_det(A);adjA=vb_adj(A);vb_need(A*adjA==detA*matid(k),"cofactor inverse identity");
 U=SAVED_CC[c]*T;expected=rho*SAVED_NN[c]*(U~*SAVED_G*U-(U~*SAVED_G*R)*adjA*(R~*SAVED_G*U)/detA);
 vb_need(den==denominator(expected)&&form==den*expected,"all quotient form entries independently reconstructed");
 vb_need(form==form~&&denominator(form)==1&&bound>0,"integer quotient and bound");1
};
vb_main()={
 my(anchor=SAVED_FIRST[1],m=SAVED_FIRST[2],rho=SAVED_FIRST[3],a=SAVED_FIRST[4],P=SAVED_FIRST[5],b2=ANCHOR_RESULT[5],num,den,
 acts2=Map(),acts3=Map(),acts4=Map(),seen2=Map(),seen3=Map(),forms2=List(),forms3=Map(),forms4=Map(),
 counts=vector(13),A,w,id,key,k,r,t,parent,c,cl,meta,stage,nq,Q,V,nv,qid,expected,Rraw,rawdelta,source,Rcan,status,stored,record,bound,
 timer=getwalltime(),fd=fileopen(Str(AUDIT_OUT,"/coverage.csv"),"w"));
 vb_need(SAVED_G==ORDINARY_G&&SAVED_G==SAVED_G~&&SAVED_W^2-SAVED_W+6*matid(22)==0,"ambient input identity");
 vb_need(SAVED_W~*SAVED_G==SAVED_G*(matid(22)-SAVED_W),"omega adjoint identity");
 vb_need(SAVED_NN==[1,2,2]&&#SAVED_CC==3,"ideal norms");
 VB_CI=vector(3,c,SAVED_CC[c]^-1);
 for(c=1,3,vb_need(VB_CI[c]*SAVED_CC[c]==matid(22)&&vb_equal_module(SAVED_CC[c],vb_twist(matid(22),c)),"whole cusp module equality"));
 vb_need(ANCHOR_RESULT[1..4]==SAVED_FIRST[1..4]&&vb_delta(P)==rho&&vb_primitive(P)&&a>=0&&a<=2,"first module metadata");
 vb_need((SAVED_X~*SAVED_G*SAVED_X)/2==m&&vb_in_module(P,Mat([SAVED_X])),"actual first ordinary vector in P");
 k=0;for(i=1,22,if(!k&&SAVED_X[i],k=sign(SAVED_X[i])));key=sum(i=1,22,(k*SAVED_X[i]+32)*2^(6*(22-i)));k=vecsearch(ORDINARY_KEYS,key);vb_need(k>0&&ORDINARY_VALUES[k]==128*m+anchor,"first ordinary vector orbit ID");
 for(i=1,#SAVED_STAB,A=SAVED_STAB[i];vb_need(denominator(A)==1&&abs(vb_det(A))==1&&A~*SAVED_G*A==SAVED_G&&A*SAVED_W==SAVED_W*A&&vb_equal_module(P,A*P),"stabilizer isometry fixing first module"));
 if(a==0,num=30613*47*rho^2;den=1029,num=46*47*rho^2;den=1);vb_need(b2^3*den<=num&&(b2+1)^3*den>num,"second exact cutoff");
 for(i=1,#AUDIT_R_IDS,id=AUDIT_R_IDS[i];w=eval(Str("R",id));key=[w[1],w[2]];vb_need(!mapisdefined(acts2,key),"rank2 disjoint kept actions");mapput(acts2,key,[4,id,w]));
 for(i=1,#AUDIT_T_IDS,id=AUDIT_T_IDS[i];w=eval(Str("T",id));key=[w[1],w[2]];vb_need(!mapisdefined(acts3,key),"rank3 disjoint kept actions");mapput(acts3,key,[4,id,w]));
 for(i=1,#AUDIT_C_IDS,id=AUDIT_C_IDS[i];w=eval(Str("C",id));key=[w[1],w[2]];vb_need(!mapisdefined(acts4,key),"rank4 disjoint actions");mapput(acts4,key,[4,id,w]));
 for(i=1,#AUDIT_W2_IDS,id=AUDIT_W2_IDS[i];w=eval(Str("W2_",id));key=[w[2],w[3]];vb_need(!mapisdefined(acts2,key),"rank2 disjoint drop actions");mapput(acts2,key,[w[1],id,w]));
 for(i=1,#AUDIT_W3_IDS,id=AUDIT_W3_IDS[i];w=eval(Str("W3_",id));key=[w[2],w[3]];vb_need(!mapisdefined(acts3,key),"rank3 disjoint drop actions");mapput(acts3,key,[w[1],id,w]));
 filewrite(fd,"rank,quotient_id,vector_index,decision,record_id");
 for(qid=1,#AUDIT_MANIFEST,
  record=AUDIT_MANIFEST[qid];vb_need(record[1]==qid,"contiguous quotient IDs");nq=record[2];bound=record[3];stage=record[5];meta=eval(Str("Q",qid));Q=eval(Str("AUDIT_Q",qid));V=eval(Str("V",qid));
  vb_need(matsize(Q)==[nq,nq],"exact input dimension");vb_quotient(qid,meta,Q,bound);counts[13]++;
  source=meta[1];rawdelta=meta[2];c=meta[3];
  if(stage==2,
   vb_need(record[6]==anchor&&record[7]==c&&record[8]==rho&&source==P&&bound==2*b2*meta[5],"second manifest and cutoff");listput(forms2,c);cl=(a+c-1)%3;
  ,if(stage==3,
   parent=record[6];w=eval(Str("R",parent));vb_need(source==w[5]&&rawdelta==w[3]&&record[7]==w[3]&&record[8]==w[4]&&record[9]==c,"third parent coverage");vb_need(bound==2*sqrtint((23*47*w[3])\5)*meta[5],"third exact cutoff");key=[parent,c];vb_need(!mapisdefined(forms3,key),"duplicate third form");mapput(forms3,key,1);cl=(w[4]+c-1)%3;
  ,vb_need(stage==4,"known stage");parent=record[6];w=eval(Str("T",parent));vb_need(source==w[5]&&rawdelta==w[3]&&record[7]==w[3]&&record[8]==w[4]&&record[9]==c,"fourth parent coverage");cl=(w[4]+c-1)%3;k=if(cl==0,47,if(m<=10,23,0));vb_need(k>0&&record[10]==k&&bound==2*k*meta[5],"fourth exact cutoff");key=[parent,c];vb_need(!mapisdefined(forms4,key),"duplicate fourth form");mapput(forms4,key,1)));
  nv=matsize(V)[2];vb_need(2*nv==record[4]&&(!nv||matsize(V)[1]==nq)&&denominator(V)==1,"saved enumeration count and dimensions");A=Map();
  for(j=1,nv,
   t=V[,j];k=0;for(i=1,#t,if(!k&&t[i],k=sign(t[i])));vb_need(k!=0,"nonzero quotient vector");key=k*t;vb_need(!mapisdefined(A,key),"unique enumeration sign pair");mapput(A,key,1);k=t~*Q*t;vb_need(k>0&&k<=bound,"exact saved vector norm");
   Rraw=concat(source,vb_line(SAVED_CC[c]*meta[4]*t,c));vb_need(denominator(Rraw)==1,"integer raw extension");expected=k/(2*meta[5]);vb_need(expected==vb_delta(Rraw),"direct raw determinant identity");
   key=[qid,j];if(stage==2,vb_need(mapisdefined(acts2,key,&stored),"rank2 full vector coverage");counts[1]++,if(stage==3,vb_need(mapisdefined(acts3,key,&stored),"rank3 full vector coverage");counts[6]++,vb_need(mapisdefined(acts4,key,&stored),"rank4 full vector coverage");counts[12]++));status=stored[1];id=stored[2];w=stored[3];
   if(stage==4,vb_need(w[3]==cl&&w[4]==Rraw,"rank4 actual candidate");filewrite(fd,Str(stage,",",qid,",",j,",4,",id));next);
   if(status==0,vb_need(w[4]==Rraw&&!vb_primitive(Rraw),"positive nonsaturation certificate");if(stage==2,counts[2]++,counts[7]++),
    vb_need(vb_primitive(Rraw),"primitive retained extension");Rcan=mathnf(Rraw);vb_need(vb_equal_module(Rraw,Rcan),"raw HNF unimodular equivalence");
    if(stage==2,
     if(status==1,vb_need(w[4]==Rcan,"rank2 ordinary module");vb_ordinary(Rcan,w[5],m,anchor);counts[3]++,
      Rcan=vb_canonical(Rcan);
      if(status==2,vb_need(vb_equal_module(w[4],Rraw)&&w[5]==Rcan&&mapisdefined(seen2,Rcan,&k)&&k==cl,"rank2 duplicate previously retained equivalent module");counts[4]++,
       vb_need(status==4&&w[3]==expected&&w[4]==cl&&w[5]==Rcan&&!mapisdefined(seen2,Rcan),"kept rank2 determinant class canonical image");mapput(seen2,Rcan,cl);counts[5]++)),
     if(status==1,vb_need(w[4]==Rcan,"rank3 ordinary module");vb_ordinary(Rcan,w[5],m,anchor);counts[8]++,
      if(status==3,vb_need(w[4]==source&&w[5]==Rcan,"second-order source module");vb_second(P,source,Rcan,rawdelta,w[6]);counts[9]++,
       if(status==2,vb_need(w[4]==Rcan&&mapisdefined(seen3,Rcan,&k)&&k==cl,"rank3 duplicate previously retained same module");counts[10]++,
        vb_need(status==4&&w[3]==expected&&w[4]==cl&&w[5]==Rcan&&!mapisdefined(seen3,Rcan),"kept rank3 determinant class exact image");mapput(seen3,Rcan,cl);counts[11]++))))
   );filewrite(fd,Str(stage,",",qid,",",j,",",status,",",id));
  );
 );
 vb_need(vecsort(Vec(forms2))==[1,2,3],"all three second classes");
 for(i=1,#AUDIT_R_IDS,for(c=1,3,vb_need(mapisdefined(forms3,[AUDIT_R_IDS[i],c]),"all third classes for every kept rank2")));
 for(i=1,#AUDIT_T_IDS,w=eval(Str("T",AUDIT_T_IDS[i]));for(c=1,3,cl=(w[4]+c-1)%3;vb_need(mapisdefined(forms4,[AUDIT_T_IDS[i],c])==(cl==0||m<=10),"all required fourth classes for every kept rank3")));
 vb_need(counts==COUNTS,"all saved COUNTS match independently verified coverage");vb_need(#acts2==counts[1]&&#acts3==counts[6]&&#acts4==counts[12],"no extraneous candidate actions");fileclose(fd);
 print("D4_INDEPENDENT_CERTIFICATE_VERIFIED anchor=",anchor," counts=",counts," ms=",getwalltime()-timer);write(Str(AUDIT_OUT,"/verification_result.gp"),"VERIFIED_COUNTS=",counts,";VERIFIED_ANCHOR=",anchor,";");1
};
vb_main();
