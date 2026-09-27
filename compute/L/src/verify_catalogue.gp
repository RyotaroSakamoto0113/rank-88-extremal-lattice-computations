default(nbthreads,1);
default(parisizemax,2000000000);
need(b,s)={if(!b,error(s))};
main()={
 my(start=getwalltime(),r,Q,sv,au,fd,U,cut);
 read(Str(ROOT,"/references/catalogue/data.gp"));fd=fileopen(Str(OUT,"/exact_forms.txt"),"w");
 for(i=1,#CATALOGUE,r=CATALOGUE[i];Q=r[2];need(Q==Q~ && matsnf(Q)==Vecrev(r[3]),"Smith divisors");sv=qfminim(Q,r[4],,2);need(sv[1]==r[5] && sv[2]==r[4],"minimum kissing number");for(j=1,#r[7],need(r[7][j]*Q*r[7][j]~==Q || r[7][j]~*Q*r[7][j]==Q,"catalogue generator isometry"));au=qfauto(Q);need(au[1]==r[6],"full group order");
  U=qflllgram(Q);need(abs(matdet(U))==1,"reduction unimodular");Q=U~*Q*U;filewrite(fd,Str(22," ",r[4]));for(j=1,22,filewrite(fd,strjoin(vector(22,k,Str(Q[j,k]))," ")));
  write(Str(OUT,"/verification.jsonl"),Str("{\"name\":\"",r[1],"\",\"minimum\":",r[4],",\"kissing_number\":",r[5],",\"group_order\":",au[1],",\"smith\":",r[3],"}"));print("CATALOGUE_VERIFIED ",r[1]);
 );fileclose(fd);print("CATALOGUE_COMPLETE milliseconds=",getwalltime()-start);
};
main();
quit;
