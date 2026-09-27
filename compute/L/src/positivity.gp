default(nbthreads,1);
default(parisizemax,4000000000);
need(b,s)={if(!b,error(s))};
main()={
 my(start=getwalltime(),eta=Mod(x,polcyclo(23)),t=eta+eta^-1,s,w,a,delta,K,J,Jb,cp,b,pol,z,cy,bs,C,H,T,B,S,U,iso,autos,G7,Barnes,Hb,found=0,bc,sc,V,A);
 read(Str(ROOT,"/workspace/rank88_exact_data.gp"));K=nfinit(polcyclo(23));s=sum(i=1,22,kronecker(i,23)*eta^i);w=(1+s)/2;a=A23(t);delta=(2-t)^5;
 need(polsturm(minpoly(a),-oo,0)==0 && polsturm(minpoly(delta),-oo,0)==0,"exact total positivity by Sturm");
 J=idealhnf(K,47,eta-21);Jb=idealhnf(K,47,eta^-1-21);need(idealhnf(K,a)==idealmul(K,J,Jb),"a generates J Jbar");
 need(idealhnf(K,delta)==idealdiv(K,K.diff,idealhnf(K,s)),"relative different23");
 cp=polcompositum(x^2-x+6,polcyclo(5),1)[1];K=nfinit(cp[1]);w=cp[2];z=cp[3];s=2*w-1;delta=(11*s*(z-z^-1)*(z+z^-1)-2*s*(z-z^-1)-115*(z+z^-1)+115)/2;
 need(polsturm(minpoly(delta),-oo,0)==0 && polsturm(minpoly(delta))==poldegree(minpoly(delta)),"exact delta4 positivity");
 write(Str(OUT,"/positivity.gp"),Str("DELTA4_MINPOL=",minpoly(delta),";A_MINPOL=",minpoly(a),";DELTA11_MINPOL=",minpoly((2-t)^5),";"));print("POSITIVITY_AND_IDEALS_VERIFIED");
 print("POSITIVITY_COMPLETE milliseconds=",getwalltime()-start);
};
main();
quit;
