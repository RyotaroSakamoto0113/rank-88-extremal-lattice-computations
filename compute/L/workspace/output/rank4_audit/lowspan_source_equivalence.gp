default(parisize,400000000);
w=Mod(x,x^2-x+6);d=2*w-1;
read("output/rank4_audit/independent_charts.gp");
read("output/rank4_audit/lowspan_data.gp");
barf(t)=subst(lift(t),x,1-w);
keymat(H)=vector(32,i,polcoef(lift(H[(i-1)%4+1,((i-1)\4)%4+1]),(i-1)\16));
assert(b,s)={if(!b,error(s))};
main()={
my(start=getabstime(),base=HH[1],done=vector(#HH),out=fileopen("output/rank4_audit/lowspan_source_equivalence_data.gp","w"));
for(j=1,#HH,for(i=1,#aut,my(A=aut[i]^-1);if(A*base*barf(A~)==HH[j],assert(denominator(lift(A))==1 && matdet(A)*barf(matdet(A))==1,"integral source equivalence");filewrite(out,Str("maps[",j,"]=",A,";"));done[j]=i;break)));
fileclose(out);
print("base_id=",ids[1]," target_ids=",ids," source_automorphism_indices=",done);
print("all_equivalent=",vecmin(done)>0," wall_ms=",getabstime()-start);
};
main();quit;
