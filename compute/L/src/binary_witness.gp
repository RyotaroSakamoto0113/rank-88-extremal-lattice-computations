default(nbthreads,1);
default(parisizemax,1000000000);
need(b,s)={if(!b,error(s))};
main()={
 my(start=getwalltime(),R,Q,V,A,C,D,vs,fs,e,f,H0,H1,hh,K,B,found=0);
 read(Str(ROOT,"/workspace/output/d3/benchmark/ambient.gp"));read(Str(ROOT,"/results/additional/binary_input.gp"));
 R=Mat([UV[1]~,W*UV[1]~,UV[2]~,W*UV[2]~]);Q=R~*G*R;A=matinverseimage(R,W*R);need(R*A==W*R,"binary O action");
 C=mathnf(concat(2*matid(4),A-matid(4)))/2;D=mathnf(concat(2*matid(4),A));
 vs=qfminim(C~*Q*C,8,,2)[3];fs=qfminim(D~*Q*D,240,,2)[3];H0=(11*G+W~*G)/23;H1=(G-2*W~*G)/23;
 for(i=1,matsize(vs)[2],e=R*C*vs[,i];if(e~*G*e!=8,next);
  for(j=1,matsize(fs)[2],f=R*D*fs[,j];if(f~*G*f!=240,next);
   forstep(eps=-1,1,2, f*=eps; if(e~*H0*f==-22 && e~*H1*f==2,
    K=Mat([2*e,W*e,f,(matid(22)-W)*f/2]);need(denominator(K)==1 && mathnf(K)==mathnf(R),"pseudobasis generates same S");
    write(Str(OUT,"/witness.gp"),Str("e=",e,";f=",f,";S_BASIS=",K,";"));
    write(Str(OUT,"/witness.json"),Str("{\"e\":",vector(22,k,Str(e[k])),",\"f\":",vector(22,k,Str(f[k])),",\"ideal_e\":\"(2,alpha)\",\"ideal_f\":\"(1,(1-alpha)/2)\",\"gram\":[[[4,0],[-22,2]],[[-20,-2],[120,0]]],\"determinant\":16}"));found=1;break(3)
   );f*=eps)
  )
 );need(found,"displayed binary pseudobasis found");print("BINARY_WITNESS_COMPLETE milliseconds=",getwalltime()-start);
};
main();
quit;
