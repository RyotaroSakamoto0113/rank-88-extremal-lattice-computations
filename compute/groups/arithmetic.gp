default(parisizemax,1000000000);
need(b,s)={if(!b,error(s))};
e=Mod(x,polcyclo(23));K=nfinit(polcyclo(23));a=1+sum(i=1,11,e^[1,2,3,4,6,8,9,12,13,16,18][i]);d=2*a-1;t=e+e^-1;AA=t^10+3*t^9-6*t^8-27*t^7-2*t^6+86*t^5+66*t^4-116*t^3-112*t^2+54*t+49;DD=(2-t)^5;JJ=idealhnf(K,47,e-21);JB=idealhnf(K,47,e^-1-21);need(idealhnf(K,AA)==idealmul(K,JJ,JB),"JJbar");need(idealhnf(K,DD)==idealdiv(K,K.diff,idealhnf(K,d)),"relative different");need(polsturm(minpoly(AA),-oo,0)==0 && polsturm(minpoly(AA))==11,"positive A");print("ARITHMETIC verified; minpoly A=",minpoly(AA));
quit;
