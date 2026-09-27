\\ Regenerate the complete q <= 12 shell, keeping file descriptors open.
\\ The vector order and file format agree with the original write-per-row code.
default(parisizemax,2000000000);
read("output/d3/benchmark/ambient.gp");
read("output/rank4_audit/ambient_qfauto.gp");
t0=getwalltime();S=qfminim(G,24);print("signed_count=",S[1]," pairs=",matsize(S[3])[2]," msec=",getwalltime()-t0);
fd=fileopen("output/rank4_audit/M_shell12.txt","w");
filewrite(fd,Str(matsize(S[3])[2]," 22"));
for(j=1,matsize(S[3])[2],filewrite(fd,Str(vector(22,i,S[3][i,j]))));
fileclose(fd);
fd=fileopen("output/rank4_audit/M_integer_data.txt","w");
filewrite(fd,"22");
for(i=1,22,filewrite(fd,Str(vector(22,j,G[i,j]))));
for(i=1,22,filewrite(fd,Str(vector(22,j,W[i,j]))));
filewrite(fd,Str(#AUT[2]));
for(k=1,#AUT[2],for(i=1,22,filewrite(fd,Str(vector(22,j,AUT[2][k][i,j])))));
fileclose(fd);
fd=fileopen("output/rank4_audit/M_shell12_count_input.txt","w");
filewrite(fd,"22 24");
for(i=1,22,filewrite(fd,Str(vector(22,j,G[i,j]))));
fileclose(fd);
print("total_ms=",getwalltime()-t0);quit;
