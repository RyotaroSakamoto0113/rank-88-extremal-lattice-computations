
default(nbthreads,1);
need(b,s)={if(!b,error(s))};
read(Str(ROOT,"/src/explicit_beta.gp"));
main()={
 my(K=nfinit(KPOL),eta=Mod(x,KPOL),al,pb,J,Jb,I);
 al=1+sum(k=1,22,if(kronecker(k,23)==1,eta^k,0));
 pb=idealhnf(K,2,1-al);J=idealhnf(K,47,eta-21);Jb=idealhnf(K,47,eta^-1-21);
 I=idealdiv(K,idealdiv(K,J,Jb),pb);
 need(idealhnf(K,BETA)==I,"beta ideal identity");
 need(2*BETA*subst(lift(BETA),x,eta^-1)==1,"beta exact norm");
 print("SEMILINEAR_COMPLETE beta ideal identity and 2*beta*bar(beta)=1");
};
main();
quit;
