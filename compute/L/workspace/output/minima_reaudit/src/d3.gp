\\ All projective determinant <= 26, using the safe ternary constant 46.
shorter_line(R,rho)={
 my(C,Q,U,sv,y,den);
 for(c=1,3,C=twist(R,c);Q=NN[c]*C~*G*C;U=qflllgram(Q);need(abs(matdet(U))==1,"small LLL");C=C*U;Q=U~*Q*U;den=denominator(Q);sv=qfminim(den*Q,den*(2*rho-2),,2);
  if(sv[1],y=C*sv[3][,1];need(NN[c]*(y~*G*y)/2<rho && denominator(linebasis(y,c))==1,"short line witness");return([c,y,NN[c]*(y~*G*y)/2])));
 []
};
main()={
 my(first=List(),reps,P,rho,Q2,sv2,Rraw,R,dr,cl,wh,key,Q3,sv3,Traw,dt,clt,sn,fdR,fdW,fdT,fdS,
 raw2=0,nonsat2=0,short2=0,duplicate2=0,keep2=0,raw3=0,nonsat3=0,sat3=0,free3=0,min2=10^9,min3=10^9,seen=Map(),record,found=List(),start=getwalltime());
 for(c=1,3,read(Str(BOOT,"/orbits",c,"/orbit_representatives.gp"));reps=MINIMA_ORBIT_REPS;for(j=1,#reps,listput(first,[c,reps[j][1],reps[j][2],CC[c]*reps[j][4]~])));
 fdR=fileopen(Str(OUT,"/rank2_modules.gp"),"w");fdW=fileopen(Str(OUT,"/rank2_shorter_witnesses.gp"),"w");fdT=fileopen(Str(OUT,"/ternary_candidates.gp"),"w");
 for(f=1,#first,record=first[f];rho=record[3];P=linebasis(record[4],record[1]);need(isprimitive(P) && delta(P)==rho,"first line saturation/determinant");P=mathnf(P);
  for(c2=1,3,Q2=quotient(P,rho,c2);sv2=qfminim(Q2[1],2*sqrtint((23*26*rho)\5)*Q2[3],,2);emit(Q2[1],2*sqrtint((23*26*rho)\5)*Q2[3],sv2,Str("second,",f,",",c2,",",rho));
   for(j=1,matsize(sv2[3])[2],raw2++;dr=(sv2[3][,j]~*Q2[1]*sv2[3][,j])/(2*Q2[3]);
    Rraw=concat(P,linebasis(CC[c2]*Q2[2]*sv2[3][,j],c2));need(denominator(Rraw)==1 && dr==delta(Rraw),"binary determinant identity");
    if(!isprimitive(Rraw),nonsat2++;next);R=mathnf(Rraw);min2=min(min2,dr);cl=(record[1]-1+c2-1)%3;
    if(mapisdefined(seen,R,&key),need(key==cl,"Steinitz path mismatch");duplicate2++;next);
    wh=shorter_line(R,rho);if(#wh,short2++;filewrite(fdW,Str("W",raw2,"=[",f,",",c2,",",j,",",R,",",wh,"] ;"));next);
    mapput(seen,R,cl);keep2++;filewrite(fdR,Str("R",keep2,"=[",f,",",record[1],",",record[2],",",rho,",",c2,",",dr,",",cl,",",R,"];"));
    for(c3=1,3,Q3=quotient(R,dr,c3);sv3=qfminim(Q3[1],52*Q3[3],,2);emit(Q3[1],52*Q3[3],sv3,Str("third,",keep2,",",dr,",",cl,",",c3));
     for(k=1,matsize(sv3[3])[2],raw3++;dt=(sv3[3][,k]~*Q3[1]*sv3[3][,k])/(2*Q3[3]);
      Traw=concat(R,linebasis(CC[c3]*Q3[2]*sv3[3][,k],c3));need(denominator(Traw)==1 && dt==delta(Traw),"ternary determinant identity");sn=isprimitive(Traw);clt=(cl+c3-1)%3;
      filewrite(fdT,Str("T",raw3,"=[",keep2,",",c3,",",k,",",dt,",",clt,",",sn,",",Traw,"];"));
      if(!sn,nonsat3++;next);sat3++;min3=min(min3,dt);if(clt==0,free3++);listput(found,[dt,clt]);
     )
    )
   )
  );
  print("FIRST_COMPLETE ",f,"/",#first," raw2=",raw2," kept=",keep2," final=",sat3," ms=",getwalltime()-start)
 );
 fileclose(fdR);fileclose(fdW);fileclose(fdT);closecert();
 fdS=fileopen(Str(OUT,"/result.gp"),"w");filewrite(fdS,Str("RESULT=[",#first,",",raw2,",",nonsat2,",",short2,",",duplicate2,",",keep2,",",raw3,",",nonsat3,",",sat3,",",free3,",",min2,",",min3,",",batchno,"];FOUND=",Vec(found),";FIRST=",Vec(first),";"));fileclose(fdS);
 print("D3_ALL_PROJECTIVE_COMPLETE result=",[#first,raw2,nonsat2,short2,duplicate2,keep2,raw3,nonsat3,sat3,free3,min2,min3,batchno]);
};
main();quit;
