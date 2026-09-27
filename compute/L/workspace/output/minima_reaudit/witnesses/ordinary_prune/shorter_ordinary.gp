\\ Load ordinary_lookup.gp first. All orbit IDs are zero based.
\\ A [] result is only "no pruning witness found"; no absence claim is made.

ordinary_key(v)={
  my(k=0,s=0,z);
  if(#v!=22,error("ordinary_key: dimension must be22"));
  for(i=1,22,if(type(v[i])!="t_INT",error("ordinary_key: noninteger coordinate"));if(!s&&v[i],s=sign(v[i])));
  if(!s,error("ordinary_key: zero vector"));
  for(i=1,22,z=s*v[i];if(z< -32||z>=32,error("ordinary_key: coordinate guard"));k=64*k+z+32);
  return(k);
};

ordinary_orbit(v,q)={
  my(k=ordinary_key(v),j=vecsearch(ORDINARY_KEYS,k),val);
  if(!j,error("ordinary_orbit: short vector is absent from complete reference sphere"));
  val=ORDINARY_VALUES[j];
  if(val\128!=q,error("ordinary_orbit: norm mismatch"));
  return(val%128);
};

ordinary_input_check(G)={
  if(G!=ORDINARY_G,error("ordinary_input_check: ambient Gram basis mismatch"));
  if(#ORDINARY_KEYS!=ORDINARY_PAIR_COUNT||#ORDINARY_VALUES!=ORDINARY_PAIR_COUNT,error("ordinary_input_check: lookup count mismatch"));
  return(1);
};

shorter_ordinary(R,m,anchor_id)={
  my(sz=matsize(R),Q,L,z,v,q,oid,first,checked=0);
  if(sz[1]!=22||sz[2]<1||sz[2]>8,error("shorter_ordinary: expected22 by1..8 integer basis"));
  if(type(m)!="t_INT"||m<1||m>ORDINARY_MAX_NORM,error("shorter_ordinary: norm cutoff outside lookup"));
  if(type(anchor_id)!="t_INT"||anchor_id<0||anchor_id>=ORDINARY_ORBIT_COUNT,error("shorter_ordinary: unknown anchor orbit"));
  for(i=1,22,for(j=1,sz[2],if(type(R[i,j])!="t_INT",error("shorter_ordinary: basis must be integral"))));
  Q=R~*ORDINARY_G*R;
  \\ qfminim is used only to propose witnesses. Every returned witness is
  \\ recomputed below with exact integer arithmetic in the original J basis.
  L=qfminim(Q,2*m,,2)[3];
  for(j=1,matsize(L)[2],
    z=L[,j];v=R*z;q=(v~*ORDINARY_G*v)/2;checked++;
    if(type(q)!="t_INT"||q<=0||q>m,error("shorter_ordinary: enumerator norm mismatch"));
    oid=ordinary_orbit(v,q);
    if(q<m||(q==m&&oid<anchor_id),
      first=0;for(i=1,22,if(!first&&v[i],first=sign(v[i])));
      v=first*v;z=first*z;
      ORDINARY_LAST_CHECK=[1,matsize(L)[2],checked];
      return([q,oid,v,z]);
    );
  );
  ORDINARY_LAST_CHECK=[0,matsize(L)[2],checked];
  return([]);
};
