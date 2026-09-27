\\ Probe saved ternary prefixes; witness search only, not a lower bound.
read(Str(OUT,"/ambient.gp"));read(Str(OUT,"/generators.gp"));
read(Str(SOURCE,"/structure.gp"));
read(Str(SOURCE,"/result.gp"));read(Str(SOURCE,"/rank3_modules.gp"));
main()={my(best=49,k,R,dt,cl3,Q,sv,Y,d,cl,fd,seen=Map(),num=0);
 fd=fileopen(Str(OUT,"/witnesses.gp"),"w");
 for(i=1,COUNTS[11],k=eval(Str("T",i));dt=k[3];cl3=k[4];R=k[5];
  for(c=1,3,cl=(cl3+c-1)%3;if(cl==0,next);
   Q=quotient(R,dt,c);sv=qfminim(Q[1],2*min(48,best)*Q[3],,2);
   for(j=1,matsize(sv[3])[2],Y=concat(R,linebasis(CC[c]*Q[2]*sv[3][,j],c));d=delta(Y);
    need(d==sv[3][,j]~*Q[1]*sv[3][,j]/(2*Q[3]),"witness determinant identity");
    if(d<=best,best=d;num++;filewrite(fd,Str("W",num,"=[",d,",",cl,",",isprimitive(Y),",",mathnf(Y),"];"));print("FOUND determinant=",d," class=",cl," saturated=",isprimitive(Y)," prefix=",i," cusp=",c))
   )
  )
 );fileclose(fd);closecert();print("PROBE_COMPLETE best=",best," witnesses=",num)
};
main();quit;
