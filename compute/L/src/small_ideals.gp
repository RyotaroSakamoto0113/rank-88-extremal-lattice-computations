default(nbthreads,1);
need(b,s)={if(!b,error(s))};
main()={
 my(F=bnfinit(x^2-x+6,1),ps,cls,z=Mod(x,polcyclo(7)),s,t,B,S);
 need(bnfcertify(F)==1 && F.cyc==[3] && F.tu[1]==2,"quadratic class and units");
 forprime(p=2,3,ps=idealprimedec(F,p);need(#ps==2 && ps[1].e==1 && ps[2].e==1 && ps[1].f==1 && ps[2].f==1,"split primes");cls=vector(2,i,bnfisprincipal(F,ps[i],0)[1]);need(Set(cls)==[1,2],"two nontrivial prime ideal classes");write(Str(OUT,"/certificate.gp"),Str("PRIMES",p,"=",ps,";CLASSES",p,"=",cls,";")));
 s=sum(i=1,6,kronecker(i,7)*z^i);t=z+z^-1;B=bnfinit(minpoly(t),1);need(bnfcertify(B)==1,"Barnes real units certified");S=concat(matrix(3,1,i,j,1),matrix(3,2,i,j,(1-bnfsignunit(B)[i,j])/2));need(matrank(Mod(S,2))==3,"Barnes real signature rank3");
 K=nfinit(polcyclo(7));need(idealhnf(K,2-t)==idealdiv(K,K.diff,idealhnf(K,s)),"Barnes relative different");need(polsturm(minpoly(2-t),-oo,0)==0,"Barnes total positivity");
 write(Str(OUT,"/certificate.gp"),Str("F_BNF=",F,";BARNES_REAL_BNF=",B,";BARNES_SIGNATURES=",S,";"));
 print("SMALL_IDEALS_COMPLETE nontrivial_classes_at_2_and_3, Barnes_positive_weight_and_units");
};
main();
quit;
