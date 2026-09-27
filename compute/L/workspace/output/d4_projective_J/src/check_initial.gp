need(b,s)={if(!b,error(s))};
read(Str(INPUTS,"/ambient.gp"));read(Str(INPUTS,"/FIRST4.gp"));read(Str(INPUTS,"/STAB4.gp"));read(Str(INPUTS,"/ordinary_representatives.gp"));
pol=t^2-t+6;bnf=bnfinit(pol,1);nf=bnf.nf;al=Mod(t,pol);pclass=bnfisprincipal(bnf,idealhnf(nf,2,al),0)[1];
main()={my(rec,x,P,m,rho,cl,U,C,II,calc,st,A,nmat=0,nprod=0);
 need(#FIRST4==78&&#STAB4==78&&#MINIMA_ORBIT_REPS==78,"all78 initial records");
 for(j=1,78,rec=FIRST4[j];x=MINIMA_ORBIT_REPS[j][4]~;m=rec[2];rho=rec[3];cl=rec[4];P=rec[5];
  need(rec[1]==j-1&&MINIMA_ORBIT_REPS[j][1]==j-1&&(x~*G*x)/2==m&&m<=12,"first norm and index");
  need(matsize(P)==[22,2]&&denominator(P)==1&&vecmax(abs(matsnf(P)))==1,"primitive rank2 basis");
  need(matdet(P~*G*P)==23*rho^2,"first line discriminant");
  U=Mat([x,W*x]);C=matinverseimage(U,P);need(U*C==P,"first F support");
  II=idealadd(nf,idealhnf(nf,C[1,1]+al*C[2,1]),idealhnf(nf,C[1,2]+al*C[2,2]));
  calc=lift(Mod(bnfisprincipal(bnf,II,0)[1],3)/Mod(pclass,3));need(calc==cl&&idealnorm(nf,II)*m==rho,"independent first class and norm");
  st=STAB4[j];need(#Set(st)==#st&&setsearch(Set(st),matid(22)),"stabilizer identity and distinctness");
  for(i=1,#st,A=st[i];need(A~*G*A==G&&A*W==W*A&&denominator(A)==1&&abs(matdet(A))==1,"stabilizer O-isometry");need(A*x==x||A*x==-x,"stabilizer fixes sign pair");nmat++;
   for(k=1,#st,need(setsearch(Set(st),A*st[k]),"stabilizer closed under products");nprod++)
  )
 );print("INITIAL78_VERIFIED matrices=",nmat," products=",nprod);
 write(OUT_JSON,Str("{\"status\":\"verified\",\"anchors\":78,\"stabilizer_matrices\":",nmat,",\"stabilizer_products\":",nprod,",\"first_classes_recomputed_by_ideal_arithmetic\":true}"));
};
main();quit;
