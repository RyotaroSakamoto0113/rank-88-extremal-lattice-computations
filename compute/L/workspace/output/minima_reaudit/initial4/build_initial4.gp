\\ Set AMBIENT, ORBITS and OUT before reading this file.
\\ Classes: 0=O, 1=p=(2,omega), 2=pbar=(2,omega-1).
need(b,s)={if(!b,error(s))};
read(AMBIENT); read(ORBITS);
need(G==G~ && matdet(G)==23^11,"ambient Gram");
need(W^2-W+6*matid(22)==0 && W~*G==G*(matid(22)-W),"omega action");
rows(M)=vector(matsize(M)[1],i,vector(matsize(M)[2],j,M[i,j]));
coords(P,K)={my(n=matsize(P)[1],A,C);for(i=1,n,for(j=i+1,n,A=matrix(2,2,a,b,P[if(a==1,i,j),b]);if(matdet(A),C=A^-1*matrix(2,matsize(K)[2],a,b,K[if(a==1,i,j),b]);need(P*C==K,"span coordinate mismatch");return(C))));error("rank less than two")};
inside(P,y)={denominator(coords(P,Mat(y)))==1};
minorcontent(P)={my(d=0,n=matsize(P)[1]);for(i=1,n,for(j=i+1,n,d=gcd(d,P[i,1]*P[j,2]-P[j,1]*P[i,2])));abs(d)};
FIRST4=List();
write(Str(OUT,"/ambient.json"),"[",rows(G),",",rows(W),"]");
gettime();
{
for(k=1,#MINIMA_ORBIT_REPS,
 my(rec=MINIMA_ORBIT_REPS[k],id=rec[1],q=rec[2],x=rec[4]~,K,P,C,idx,D,rho,cl=0,lam=-1,indexcl=-1,y=vector(22,i,0)~,B,vden=1,BC,WP);
 K=Mat([x,W*x]);
 need(x~*G*x==2*q,"anchor norm");
 P=mathnf(matkerint(matkerint(K~)~));
 need(matsize(P)==[22,2] && minorcontent(P)==1,"primitive saturation basis");
 C=coords(P,K); need(denominator(C)==1,"Ox inclusion");
 idx=abs(matdet(C)); need(idx==minorcontent(K) && (idx==1 || idx==2),"anchor index");
 D=matdet(P~*G*P)/23; rho=sqrtint(D);need(rho^2==D && q==idx*rho,"line determinant");
 WP=coords(P,W*P);need(denominator(WP)==1,"omega stability");
 if(idx==1,
   B=K,
   need(q==12 && rho==6,"nonprincipal anchor norm");
   for(j=1,2,if(!inside(K,P[,j]),y=P[,j];break));
   need(y!=0 && !inside(K,y) && inside(K,2*y),"quotient generator");
   lam=if(inside(K,W*y),0,1);
   need(inside(K,W*y-lam*y) && !inside(K,W*y-(1-lam)*y),"quotient omega eigenvalue");
   indexcl=if(lam==0,1,2);cl=3-indexcl;vden=2;
   B=if(cl==1,Mat([x,W*x/2]),Mat([x,(W-matid(22))*x/2]));
 );
 need(denominator(B)==1,"explicit ideal basis integrality");
 BC=coords(P,B);need(denominator(BC)==1 && abs(matdet(BC))==1,"explicit ideal basis equality");
 need(cl==0 || cl+indexcl==3,"inverse index class");
 listput(FIRST4,[id,q,rho,cl,P]);
 write(Str(OUT,"/records.jsonl"),[id,q,rho,cl,idx,lam,indexcl,Vec(x),rows(P),rows(K),rows(C),Vec(y),rows(B),vden,rows(BC),rows(WP)]);
);
}
initial4_ms=gettime();
write(Str(OUT,"/FIRST4.gp"),"FIRST4=",Vec(FIRST4),";");
write(Str(OUT,"/gp_timing.json"),"{\"gp_cpu_ms\":",initial4_ms,",\"records\":",#FIRST4,"}");
print("verified initial lines: ",#FIRST4,"; GP CPU ms: ",initial4_ms);
