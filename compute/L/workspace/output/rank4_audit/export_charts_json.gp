w=Mod(x,x^2-x+6);
read("output/rank4_audit/independent_charts.gp");
rowkey(v)={my(a=vector(8,i,polcoef(lift(v[1+(i-1)%4]),(i-1)\4)));if(lex(a,-a)<0,a,-a)};
chartkey(R)={vecsort(vector(4,i,rowkey(R[i,])))};
matkey(U)=vector(32,i,polcoef(lift(U[(i-1)%4+1,((i-1)\4)%4+1]),(i-1)\16));
KS=Set();RS=List();
for(i=1,#RR,kk=chartkey(RR[i]);if(!setsearch(KS,kk),KS=setunion(KS,Set([kk]));listput(RS,matkey(23*RR[i]))));
fd=fileopen("output/rank4_audit/charts_scaled23.json","w");filewrite(fd,Vec(RS));fileclose(fd);quit;
