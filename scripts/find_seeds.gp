\\ Exploratory fixed-seed LLL; only exact norm checks enter the certificate.
need(b,s)=if(!b,error(s));
main()={
 my(n=88,U0,H0,R,U,H,perm,v,q,fd,ft,seen=Map(),count=0,added,start,cpu,tlll,clll,lllw,lllc,preCpu,ii,jj,cc,key);
 need(matsize(G)==[n,n] && denominator(G)==1 && G==G~,"integer Gram input");
 setrand(RNGSEED);fd=fileopen(Str(OUT,"/candidates.jsonl"),"w");ft=fileopen(Str(OUT,"/trials.csv"),"w");filewrite(ft,"trial,wall_ms,cpu_ms,lll_wall_ms,lll_cpu_ms,minimum_basis_norm,new_candidates,total_candidates");
 start=getwalltime();gettime();U0=qflllgram(G);H0=U0~*G*U0;need(abs(matdet(U0))==1,"initial LLL basis");
 print("INITIAL_LLL ",getwalltime()-start," ms min=",vecmin(vector(n,i,H0[i,i])));
 for(trial=0,TRIALS,
  start=getwalltime();gettime();
  if(trial==0,U=U0;H=H0;lllw=0;lllc=0;preCpu=0,
   perm=vector(n,i,i);forstep(i=n,2,-1,jj=1+random(i);cc=perm[i];perm[i]=perm[jj];perm[jj]=cc);
   R=matrix(n,n,i,j,i==perm[j]);
   for(s=1,n,ii=1+random(n);jj=1+random(n-1);if(jj>=ii,jj++);cc=if(random(2),1,-1);R[,jj]+=cc*R[,ii]);
   H=R~*H0*R;tlll=getwalltime();preCpu=gettime();U=qflllgram(H);lllw=getwalltime()-tlll;lllc=gettime();U=U0*R*U;H=U~*G*U;
   need(abs(matdet(U))==1,"randomized LLL basis")
  );
  added=0;
  for(i=1,n,
   q=H[i,i];need(q>=8,"found norm below assumed minimum8");
   if(q==8,v=U[,i];for(k=1,n,if(v[k],if(v[k]<0,v=-v);break));if(!mapisdefined(seen,v),mapput(seen,v,1);need(v~*G*v==8,"exact seed norm");filewrite(fd,Str("{\"trial\":",trial,",\"kind\":\"basis\",\"x\":",Vec(v),"}"));added++;count++));
   for(j=i+1,n,forstep(sgn=-1,1,2,
    if(H[i,i]+H[j,j]+2*sgn*H[i,j]==8,v=U[,i]+sgn*U[,j];for(k=1,n,if(v[k],if(v[k]<0,v=-v);break));if(!mapisdefined(seen,v),mapput(seen,v,1);need(v~*G*v==8,"exact pair seed norm");filewrite(fd,Str("{\"trial\":",trial,",\"kind\":\"pair\",\"x\":",Vec(v),"}"));added++;count++))
   ))
  );
  cpu=gettime()+lllc+preCpu;filewrite(ft,Str(trial,",",getwalltime()-start,",",cpu,",",lllw,",",lllc,",",vecmin(vector(n,i,H[i,i])),",",added,",",count));
  print("TRIAL ",trial," candidates=",count," new=",added," min=",vecmin(vector(n,i,H[i,i]))," ms=",getwalltime()-start)
 );
 fileclose(fd);fileclose(ft);print("SEED_SEARCH_COMPLETE ",count);
};
main();quit;
