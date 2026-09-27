\\ One ordinary first-vector orbit per job. Dependencies are certified separately.
canonical_binary(R,st)={my(best=R,A);for(i=1,#st,A=mathnf(st[i]*R);if(cmp(A,best)<0,best=A));best};
second_order_witness(P,rho,R,dr,T,st)={
 my(C,K,ss,U,V,Q,den,ev,y,value,Rp,canp);
 for(c=1,3,
  C=twist(T,c);K=matinverseimage(C,twist(P,c));need(denominator(K)==1 && isprimitive(K),"small quotient kernel");
  ss=matsnf(K,1);U=ss[1]^-1;V=U[,1..4];need(abs(matdet(concat(K,V)))==1,"small quotient complement");
  C=C*V;Q=rho*NN[c]*C~*schur(P)*C;U=qflllgram(Q);need(abs(matdet(U))==1,"small quotient LLL");C=C*U;Q=U~*Q*U;den=denominator(Q);ev=qfminim(den*Q,2*dr*den,,2);
  for(j=1,matsize(ev[3])[2],y=C*ev[3][,j];value=NN[c]*rho*(y~*schur(P)*y)/2;
   need(denominator(value)==1 && value>0 && value<=dr,"second witness determinant bound");
   Rp=concat(P,linebasis(y,c));need(denominator(Rp)==1 && delta(Rp)==value,"second witness determinant identity");
   if(value<dr,return([c,value,y,Rp,0]));
   if(isprimitive(Rp),canp=canonical_binary(mathnf(Rp),st);if(cmp(canp,R)<0,return([c,value,y,Rp,canp])))
  )
 );[]
};
main()={
 my(rec=FIRST4[ANCHOR+1],m,rho,a,P,st,b2=0,numer,denom,Q2,sv2,dr,Rraw,R,cl2,ow2,canR,key,Q3,sv3,dt,Traw,T,cl3,ow3,sw,Q4,sv4,goal,cl4,
  seen2=Map(),seen3=Map(),fd2,fd3,fd4,fdw,fdmeta,fds,raw2=0,nonsat2=0,ordinary2=0,duplicate2=0,keep2=0,raw3=0,nonsat3=0,ordinary3=0,second_order3=0,duplicate3=0,keep3=0,raw4=0,start=getwalltime(),i2,i3,i4,meta2,meta3,meta4);
 need(rec[1]==ANCHOR,"first anchor ID");m=rec[2];rho=rec[3];a=rec[4];P=rec[5];st=STAB4[ANCHOR+1];
 need(m<=12 && rho>=6 && isprimitive(P) && delta(P)==rho,"first module");ordinary_input_check(G);
 fds=fileopen(Str(OUT,"/structure.gp"),"w");filewrite(fds,Str("SAVED_G=",G,";SAVED_W=",W,";SAVED_CC=",CC,";SAVED_NN=",NN,";SAVED_FIRST=",rec,";SAVED_STAB=",st,";SAVED_X=",MINIMA_ORBIT_REPS[ANCHOR+1][4]~,";"));fileclose(fds);
 numer=46*47*rho^2;denom=1;while((b2+1)^3*denom<=numer,b2++);
 fd2=fileopen(Str(OUT,"/rank2_modules.gp"),"w");fd3=fileopen(Str(OUT,"/rank3_modules.gp"),"w");fd4=fileopen(Str(OUT,"/rank4_candidates.gp"),"w");fdw=fileopen(Str(OUT,"/pruning_witnesses.gp"),"w");fdmeta=fileopen(Str(OUT,"/quotient_metadata.gp"),"w");
 for(c2=1,3,Q2=quotient(P,rho,c2);sv2=qfminim(Q2[1],2*b2*Q2[3],,2);meta2=emit(Q2[1],2*b2*Q2[3],sv2,Str("second,",ANCHOR,",",c2,",",rho));filewrite(fdmeta,Str("Q",meta2,"=[",P,",",rho,",",c2,",",Q2[2],",",Q2[3],"];"));
  for(j=1,matsize(sv2[3])[2],raw2++;dr=(sv2[3][,j]~*Q2[1]*sv2[3][,j])/(2*Q2[3]);Rraw=concat(P,linebasis(CC[c2]*Q2[2]*sv2[3][,j],c2));need(denominator(Rraw)==1 && dr==delta(Rraw),"second determinant identity");
   if(!isprimitive(Rraw),nonsat2++;filewrite(fdw,Str("W2_",raw2,"=[0,",meta2,",",j,",",Rraw,"];"));next);
   need(dr>=16,"certified d2 lower bound violated");R=mathnf(Rraw);cl2=(a+c2-1)%3;
   ow2=shorter_ordinary(R,m,ANCHOR);if(#ow2,ordinary2++;filewrite(fdw,Str("W2_",raw2,"=[1,",meta2,",",j,",",R,",",ow2,"];"));next);
   canR=canonical_binary(R,st);if(mapisdefined(seen2,canR,&key),need(key==cl2,"rank2 class inconsistency");duplicate2++;filewrite(fdw,Str("W2_",raw2,"=[2,",meta2,",",j,",",R,",",canR,"];"));next);
   R=canR;mapput(seen2,R,cl2);keep2++;i2=keep2;filewrite(fd2,Str("R",i2,"=[",meta2,",",j,",",dr,",",cl2,",",R,"];"));
   for(c3=1,3,Q3=quotient(R,dr,c3);sv3=qfminim(Q3[1],2*sqrtint((23*47*dr)\5)*Q3[3],,2);meta3=emit(Q3[1],2*sqrtint((23*47*dr)\5)*Q3[3],sv3,Str("third,",i2,",",dr,",",cl2,",",c3));filewrite(fdmeta,Str("Q",meta3,"=[",R,",",dr,",",c3,",",Q3[2],",",Q3[3],"];"));
    for(k=1,matsize(sv3[3])[2],raw3++;dt=(sv3[3][,k]~*Q3[1]*sv3[3][,k])/(2*Q3[3]);Traw=concat(R,linebasis(CC[c3]*Q3[2]*sv3[3][,k],c3));need(denominator(Traw)==1 && dt==delta(Traw),"third determinant identity");
     if(!isprimitive(Traw),nonsat3++;filewrite(fdw,Str("W3_",raw3,"=[0,",meta3,",",k,",",Traw,"];"));next);
     need(dt>=27,"certified projective d3 lower bound violated");T=mathnf(Traw);cl3=(cl2+c3-1)%3;
     ow3=shorter_ordinary(T,m,ANCHOR);if(#ow3,ordinary3++;filewrite(fdw,Str("W3_",raw3,"=[1,",meta3,",",k,",",T,",",ow3,"];"));next);
     sw=second_order_witness(P,rho,R,dr,T,st);if(#sw,second_order3++;filewrite(fdw,Str("W3_",raw3,"=[3,",meta3,",",k,",",R,",",T,",",sw,"];"));next);
     if(mapisdefined(seen3,T,&key),need(key==cl3,"rank3 class inconsistency");duplicate3++;filewrite(fdw,Str("W3_",raw3,"=[2,",meta3,",",k,",",T,"];"));next);
     mapput(seen3,T,cl3);keep3++;i3=keep3;filewrite(fd3,Str("T",i3,"=[",meta3,",",k,",",dt,",",cl3,",",T,"];"));
     for(c4=1,3,cl4=(cl3+c4-1)%3;goal=47;if(!goal,next);
      Q4=quotient(T,dt,c4);sv4=qfminim(Q4[1],2*goal*Q4[3],,2);meta4=emit(Q4[1],2*goal*Q4[3],sv4,Str("fourth,",i3,",",dt,",",cl3,",",c4,",",goal));filewrite(fdmeta,Str("Q",meta4,"=[",T,",",dt,",",c4,",",Q4[2],",",Q4[3],"];"));
      for(l=1,matsize(sv4[3])[2],raw4++;filewrite(fd4,Str("C",raw4,"=[",meta4,",",l,",",cl4,",",concat(T,linebasis(CC[c4]*Q4[2]*sv4[3][,l],c4)),"];")))
     )
    )
   );
   if(keep2%50==0,print("PREFIX ",ANCHOR," R2=",keep2," raw3=",raw3," kept3=",keep3," ms=",getwalltime()-start))
  )
 );
 fileclose(fd2);fileclose(fd3);fileclose(fd4);fileclose(fdw);fileclose(fdmeta);closecert();fds=fileopen(Str(OUT,"/result.gp"),"w");filewrite(fds,Str("COUNTS=",[raw2,nonsat2,ordinary2,duplicate2,keep2,raw3,nonsat3,ordinary3,second_order3,duplicate3,keep3,raw4,batchno],";ANCHOR_RESULT=",[ANCHOR,m,rho,a,b2],";"));fileclose(fds);
 print("D4_ANCHOR_COMPLETE ",ANCHOR," counts=",[raw2,nonsat2,ordinary2,duplicate2,keep2,raw3,nonsat3,ordinary3,second_order3,duplicate3,keep3,raw4,batchno]," ms=",getwalltime()-start);
};
main();quit;
