\\ Independent full quotient audit: no shorter-line pruning, no root quotient function.
require(b,s)={if(!b,error(s))};
read(Str(BOOT,"/ambient.gp"));read(Str(BOOT,"/cusps.gp"));
N0=22;I0=matid(N0);weights=[1,2,2];
require(G==G~ && denominator(G)==1 && qfsign(G)==[22,0],"positive integral G");
require(denominator(W)==1 && W^2-W+6*I0==0 && W~*G==G*(I0-W),"O action");
require(denominator((G-2*W~*G)/23)==1 && denominator((11*G+W~*G)/23)==1,"Hermitian integral");
Cbasic=[I0,mathnf(concat(2*I0,W-I0))/2,mathnf(concat(2*I0,W))/2];
for(aa=1,3,change=Cbasic[aa]^-1*CC[aa];require(denominator(change)==1 && abs(matdet(change))==1,"cusp same lattice");require(CI[aa]*CC[aa]==I0 && GG[aa]==weights[aa]*CC[aa]~*G*CC[aa],"cusp metric"));
ind_line(y,c)={if(c==1,Mat([y,W*y]),if(c==2,Mat([2*y,W*y]),Mat([2*y,(W-I0)*y])))};
ind_sat(X)={my(sn=matsnf(X),r=matsize(X)[2]);require(denominator(X)==1,"integral submodule");sum(i=1,#sn,sn[i]!=0)==r && vecmax(abs(sn))==1};
ind_delta(X)={my(r=matsize(X)[2]/2,z=matdet(X~*G*X)/23^r,s);require(denominator(z)==1 && z>0,"trace discriminant");s=sqrtint(z);require(s*s==z,"Hermitian determinant integer");s};
ind_quot(X,detX,c)={
 my(S=G-G*X*(X~*G*X)^-1*X~*G,A=detX*weights[c]*CC[c]~*S*CC[c],dd,kt,Q,T,K,r=matsize(X)[2]);
 dd=denominator(A);kt=qflllgram(dd*A,4);K=kt[1];T=kt[2];
 require(matsize(K)==[22,r] && matsize(T)==[22,22-r],"quotient dimensions");
 require(A*K==0 && denominator(concat(K,T))==1 && abs(matdet(concat(K,T)))==1,"independent kernel-image unimodular split");
 Q=T~*A*T;dd=denominator(Q);Q*=dd;require(Q==Q~ && denominator(Q)==1 && qfsign(Q)==[22-r,0],"positive exact quotient");
 [Q,T,dd]
};
orient(v)={for(i=1,#v,if(v[i],return(if(v[i]<0,-v,v))));v};
fd=fileopen(Str(OUT,"/exact_inputs.txt"),"w");fm=fileopen(Str(OUT,"/exact_manifest.csv"),"w");fr=fileopen(Str(OUT,"/all_binary_modules.gp"),"w");ft=fileopen(Str(OUT,"/all_ternary_candidates.gp"),"w");fv=fileopen(Str(OUT,"/enumerated_vectors.gp"),"w");batch=0;
ind_emit(Q,b,sv,meta)={my(nn=matsize(Q)[1],V=sv[3],ss);require(sv[1]==2*matsize(V)[2],"signed count");ss=Set(vector(matsize(V)[2],i,orient(V[,i])));require(#ss==matsize(V)[2],"unique sign pairs");for(i=1,matsize(V)[2],require(V[,i]~*Q*V[,i]>0 && V[,i]~*Q*V[,i]<=b,"vector within bound"));batch++;filewrite(fd,Str(nn," ",b));for(i=1,nn,filewrite(fd,strjoin(vector(nn,j,Str(Q[i,j]))," ")));filewrite(fm,Str(batch,",",nn,",",b,",",sv[1],",",meta));filewrite(fv,Str("V",batch,"=",V,";"));batch};
main()={
 my(P,entry,first=List(),rho,q2,sv2,Y,dr,R,q3,sv3,Z,dt,cl,cl3,seen=Map(),raw2=0,non2=0,unique2=0,dup2=0,raw3=0,non3=0,sat3=0,free3=0,min2=10^9,min3=10^9,clk=getwalltime(),bd);
 for(c=1,3,read(Str(BOOT,"/orbits",c,"/orbit_representatives.gp"));for(i=1,#MINIMA_ORBIT_REPS,entry=MINIMA_ORBIT_REPS[i];rho=entry[2];P=ind_line(CC[c]*entry[4]~,c);require(rho<=10 && rho>=6 && ind_sat(P) && ind_delta(P)==rho,"initial primitive line");listput(first,[c,rho,mathnf(P)])));
 for(f=1,#first,P=first[f][3];rho=first[f][2];
  for(c=1,3,q2=ind_quot(P,rho,c);bd=2*sqrtint((23*26*rho)\5)*q2[3];sv2=qfminim(q2[1],bd,,2);ind_emit(q2[1],bd,sv2,Str("second,",f,",",c));
   for(j=1,matsize(sv2[3])[2],raw2++;Y=concat(P,ind_line(CC[c]*q2[2]*sv2[3][,j],c));dr=ind_delta(Y);require(dr==sv2[3][,j]~*q2[1]*sv2[3][,j]/(2*q2[3]),"binary Schur determinant");
    if(!ind_sat(Y),non2++;next);min2=min(min2,dr);R=mathnf(Y);cl=(first[f][1]+c-2)%3;
    if(mapisdefined(seen,R),require(mapget(seen,R)==cl,"class consistency");dup2++;next);mapput(seen,R,cl);unique2++;filewrite(fr,Str("R",unique2,"=[",f,",",c,",",j,",",dr,",",cl,",",R,"];"));
    for(c3=1,3,q3=ind_quot(R,dr,c3);sv3=qfminim(q3[1],52*q3[3],,2);ind_emit(q3[1],52*q3[3],sv3,Str("third,",unique2,",",c3));
     for(k=1,matsize(sv3[3])[2],raw3++;Z=concat(R,ind_line(CC[c3]*q3[2]*sv3[3][,k],c3));dt=ind_delta(Z);require(dt==sv3[3][,k]~*q3[1]*sv3[3][,k]/(2*q3[3]),"ternary Schur determinant");cl3=(cl+c3-1)%3;filewrite(ft,Str("T",raw3,"=[",unique2,",",c3,",",dt,",",cl3,",",ind_sat(Z),",",Z,"];"));if(!ind_sat(Z),non3++;next);sat3++;min3=min(min3,dt);if(cl3==0,free3++))
    )
   )
  );print("FIRST ",f,"/",#first," raw2=",raw2," unique2=",unique2," raw3=",raw3," milliseconds=",getwalltime()-clk)
 );
 fileclose(fd);fileclose(fm);fileclose(fr);fileclose(ft);fileclose(fv);
 fd=fileopen(Str(OUT,"/independent_result.json"),"w");filewrite(fd,Str("{\"first_orbits\":",#first,",\"raw_binary\":",raw2,",\"nonsaturated_binary\":",non2,",\"unique_binary\":",unique2,",\"duplicates\":",dup2,",\"raw_ternary\":",raw3,",\"nonsaturated_ternary\":",non3,",\"saturated_ternary\":",sat3,",\"free_ternary\":",free3,",\"min_binary\":",min2,",\"min_ternary\":",min3,",\"forms\":",batch,"}"));fileclose(fd);
 print("INDEPENDENT_D3_ENUMERATION_COMPLETE ",batch," forms, ",raw3," ternary candidates");
};
main();quit;
