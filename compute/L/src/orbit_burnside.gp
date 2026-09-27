default(nbthreads,1);
default(parisizemax,4000000000);
need(b,s)={if(!b,error(s))};
main()={
 my(start=getwalltime(),base=Str(ROOT,"/workspace/output/minima_reaudit"),p,auto,forms,rows,idx=0,A,V,Q,co,sgn);
 p=readstr(Str(base,"/LATEST_BOOTSTRAP.txt"))[1];read(Str(p,"/bootstrap/ambient.gp"));read(Str(p,"/bootstrap/cusps.gp"));read(Str(p,"/bootstrap/generators.gp"));
 auto=qfauto(G);need(auto[1]==12144,"full automorphism order");for(k=1,#auto[2],need(auto[2][k]*W==W*auto[2][k],"all generators O linear"));
 need(matsnf(G)==concat(vector(11,i,23),vector(11,i,1)),"Smith invariants");write(Str(OUT,"/full_automorphisms.gp"),Str("FULL_AUT=",auto,";"));print("FULL_AUT_VERIFIED order12144 O_linear Smith1^11*23^11");
 read(Str(base,"/witnesses/orbits/group_matrices.gp"));read(Str(ROOT,"/results/group/conjugacy.gp"));need(#MINIMA_GROUP==12144 && vecsum(CLASS_SIZES)==12144,"conjugacy data");
 forms=fileopen(Str(OUT,"/exact_forms.txt"),"w");rows=fileopen(Str(OUT,"/burnside_rows.jsonl"),"w");
 for(c=1,3,
  for(j=1,#CLASS_REPS,
   A=CI[c]*MINIMA_GROUP[CLASS_REPS[j]]*CC[c];need(denominator(A)==1 && A~*GG[c]*A==GG[c],"cusp isometry");
   V=matkerint(A-matid(22));need((A-matid(22))*V==0,"fixed lattice kernel");
   if(matsize(V)[2],Q=V~*GG[c]*V;U=qflllgram(Q);need(abs(matdet(U))==1,"fixed space unimodular reduction");V=V*U;Q=U~*Q*U;
    idx++;filewrite(forms,Str(matsize(Q)[1]," 42"));for(k=1,matsize(Q)[1],filewrite(forms,strjoin(vector(matsize(Q)[1],l,Str(Q[k,l]))," ")));
    co=2*Vec(qfrep(Q,21,1));,
    co=vector(21,i,0)
   );
   filewrite(rows,Str("{\"cusp\":",c,",\"class\":",j,",\"class_size\":",CLASS_SIZES[j],",\"dimension\":",matsize(V)[2],",\"form_index\":",if(matsize(V)[2],idx,0),",\"counts\":",co,"}"));
   write(Str(OUT,"/fixed_lattices.gp"),Str("FIXED_",c,"_",j,"=",V,";"));
   print("FIXED_COUNT cusp=",c," class=",j," dimension=",matsize(V)[2]," ms=",getwalltime()-start);
  );
 );fileclose(forms);fileclose(rows);print("ORBIT_BURNSIDE_COMPLETE forms=",idx," milliseconds=",getwalltime()-start);
};
main();
quit;
