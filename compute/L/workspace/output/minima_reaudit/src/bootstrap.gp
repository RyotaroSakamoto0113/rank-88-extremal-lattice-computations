\\ OUT and common.gp are supplied by the timed driver.
main()={
 my(sv,fd,fg,A,cut=20);
 for(c=1,3,
  sv=qfminim(GG[c],cut,,2);emit(GG[c],cut,sv,Str("initial_cusp,",c));
  need(sum(i=1,matsize(sv[3])[2],sv[3][,i]~*GG[c]*sv[3][,i]<12)==0,"projective norm below 6 found");
  fd=fileopen(Str(OUT,"/cusp",c,"_shell.txt"),"w");filewrite(fd,Str(matsize(sv[3])[2]," 22"));
  for(j=1,matsize(sv[3])[2],filewrite(fd,strjoin(vector(22,i,Str(sv[3][i,j]))," ")));fileclose(fd);
  fg=fileopen(Str(OUT,"/cusp",c,"_integer_data.txt"),"w");filewrite(fg,"22");
  for(i=1,22,filewrite(fg,strjoin(vector(22,j,Str(GG[c][i,j]))," ")));
  A=CI[c]*W*CC[c];need(denominator(A)==1,"cusp omega matrix");for(i=1,22,filewrite(fg,strjoin(vector(22,j,Str(A[i,j]))," ")));
  filewrite(fg,Str(#AUT[2]));for(k=1,#AUT[2],A=CI[c]*AUT[2][k]*CC[c];need(denominator(A)==1,"cusp automorphism matrix");for(i=1,22,filewrite(fg,strjoin(vector(22,j,Str(A[i,j]))," "))));fileclose(fg);
  print("CUSP_COMPLETE ",c," signed_vectors=",sv[1]);
 );
 fd=fileopen(Str(OUT,"/cusps.gp"),"w");filewrite(fd,Str("CC=",CC,";GG=",GG,";CI=",CI,";"));fileclose(fd);
 closecert();print("BOOTSTRAP_COMPLETE batches=",batchno);
};
main();quit;
