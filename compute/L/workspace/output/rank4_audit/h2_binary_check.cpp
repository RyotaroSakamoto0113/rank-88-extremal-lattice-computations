#include <gmpxx.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>
using Q=mpq_class; using Z=mpz_class;
using Clock=std::chrono::steady_clock;
struct Enumerator {
 int n=0; Q bound; std::vector<std::vector<Q>> a,u; std::vector<Q>d;
 std::vector<std::vector<long>> vectors; std::vector<long>x,first; unsigned long long nodes=0,pairs=0; bool complete=true;
 double seconds_limit; Clock::time_point started;
 Enumerator(std::istream&in,double limit):seconds_limit(limit){
  in>>n>>bound; if(!in||n<1||n>100||bound<0)throw std::runtime_error("invalid header");bound.canonicalize();
  a.assign(n,std::vector<Q>(n));u=a;d.resize(n);x.resize(n);
  for(auto&row:a)for(auto&v:row){in>>v;v.canonicalize();}
  if(!in)throw std::runtime_error("incomplete input");
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[i][j]!=a[j][i])throw std::runtime_error("not symmetric");
  for(int i=0;i<n;i++){
   d[i]=a[i][i];for(int k=0;k<i;k++)d[i]-=d[k]*u[k][i]*u[k][i];
   if(d[i]<=0)throw std::runtime_error("not positive definite");u[i][i]=1;
   for(int j=i+1;j<n;j++){Q v=a[i][j];for(int k=0;k<i;k++)v-=d[k]*u[k][i]*u[k][j];u[i][j]=v/d[i];}
  }
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){Q v=0;for(int k=0;k<n;k++)v+=d[k]*u[k][i]*u[k][j];if(v!=a[i][j])throw std::runtime_error("LDL identity failed");}
 }
 void rec(int i,const Q&remaining,bool nonzero){
  if((nodes&8191)==0 && std::chrono::duration<double>(Clock::now()-started).count()>seconds_limit){complete=false;throw std::runtime_error("time limit");}
  if(i<0){if(nonzero){pairs++; vectors.push_back(x); if(first.empty())first=x;}return;}
  Q c=0;for(int j=i+1;j<n;j++)c+=u[i][j]*x[j];
  Q radius=remaining/d[i]; Z aa=c.get_num(),bb=c.get_den();
  Z scaled=radius.get_num()*bb*bb, denom=radius.get_den(),integer_radius;
  mpz_fdiv_q(integer_radius.get_mpz_t(),scaled.get_mpz_t(),denom.get_mpz_t());
  Z k;mpz_sqrt(k.get_mpz_t(),integer_radius.get_mpz_t());
  Z lo,hi,numerator=-k-aa;mpz_cdiv_q(lo.get_mpz_t(),numerator.get_mpz_t(),bb.get_mpz_t());
  numerator=k-aa;mpz_fdiv_q(hi.get_mpz_t(),numerator.get_mpz_t(),bb.get_mpz_t());
  if(!nonzero&&lo<0)lo=0; // one representative per pair {v,-v}
  if(!lo.fits_slong_p()||!hi.fits_slong_p())throw std::runtime_error("coordinate overflow");
  const long low=lo.get_si(),high=hi.get_si();
  for(long t=low;t<=high;t++){
   x[i]=t;nodes++;Q v=c+t,rem=remaining-d[i]*v*v;
   if(rem<0)throw std::runtime_error("exact radius logic failed");
   rec(i-1,rem,nonzero||t!=0);
   if(t==high)break;
  }
 }
};

#include <array>
#include <limits>
using Int=long long;
Int narrow(__int128 v){if(v>1000000000 || v< -1000000000)throw std::runtime_error("safe integer bound exceeded");return (Int)v;}
Int dot(const std::vector<Int>&a,const std::vector<long>&b){__int128 s=0;for(size_t i=0;i<a.size();i++)s+=(__int128)a[i]*b[i];return narrow(s);}
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("usage: h2_binary_check input");
 std::ifstream in(argv[1]);if(!in)throw std::runtime_error("input missing");
 auto started=Clock::now();Enumerator e(in,120);int n=e.n;
 std::vector<std::vector<Int>> G(n,std::vector<Int>(n)),W=G,GW=G;
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){
   if(e.a[i][j].get_den()!=1 || !e.a[i][j].get_num().fits_slong_p())throw std::runtime_error("nonintegral G");
   G[i][j]=narrow(e.a[i][j].get_num().get_si());in>>W[i][j];W[i][j]=narrow(W[i][j]);
 }
 if(!in)throw std::runtime_error("bad W");
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){__int128 s=0;for(int k=0;k<n;k++)s+=(__int128)G[i][k]*W[k][j];GW[i][j]=narrow(s);}
 e.started=Clock::now();e.rec(n-1,e.bound,false);
 std::vector<std::vector<int>> shells(10);
 std::vector<std::vector<Int>> g(e.vectors.size(),std::vector<Int>(n)),gw=g;
 for(size_t k=0;k<e.vectors.size();k++){
   for(long v:e.vectors[k])narrow(v);
   for(int j=0;j<n;j++){__int128 a=0,b=0;for(int i=0;i<n;i++){a+=(__int128)e.vectors[k][i]*G[i][j];b+=(__int128)e.vectors[k][i]*GW[i][j];}g[k][j]=narrow(a);gw[k][j]=narrow(b);}
   Int q2=dot(g[k],e.vectors[k]);
   if(q2<=0 || q2>18 || q2%2)throw std::runtime_error("exact norm failure");
   if(dot(gw[k],e.vectors[k])!=q2/2)throw std::runtime_error("diagonal omega convention failure");
   shells[q2/2].push_back(k);
 }
 auto enum_done=Clock::now();
 // Only one representative of each sign pair is enumerated. The two signs
 // of y supply all ordered pairs once x has its chosen representative.
 const std::array<std::array<int,4>,4> types={{{6,9,2,2},{6,9,4,-2},{8,8,3,2},{8,8,5,-2}}};
 std::array<unsigned long long,4> counts={0,0,0,0};
 std::array<std::array<int,3>,4> witnesses;for(auto &w:witnesses)w={-1,-1,0};
 unsigned long long comparisons=0,positive_controls=0;
 // Traverse each shell pair once and check both conjugate target values.
 for(int c=0;c<2;c++){
   int a=c?8:6,b=c?8:9;
   for(int i:shells[a])for(int j:shells[b]){
     Int t=dot(g[i],e.vectors[j]),s=dot(gw[i],e.vectors[j]);comparisons++;
     // A positive control verifies exactly every h(x,x)=q(x).
     if(i==j){if(t!=2*a||s!=a)throw std::runtime_error("positive control failed");positive_controls++;}
     for(int k=2*c;k<2*c+2;k++)for(int sign:{-1,1}){
       if(sign*t==2*types[k][2]+types[k][3] && sign*s==types[k][2]+12*types[k][3]){if(counts[k]==0)witnesses[k]={i,j,sign};counts[k]++;}
     }
   }
 }
 auto ended=Clock::now();
 std::cout<<"{\"input\":\""<<argv[1]<<"\",\"complete\":true,\"dimension\":"<<n<<",\"bound\":18,\"pairs\":"<<e.vectors.size()<<",\"shell_pairs\":{";
 for(int q=1;q<=9;q++){if(q>1)std::cout<<",";std::cout<<"\""<<q<<"\":"<<shells[q].size();}
 std::cout<<"},\"ordered_sign_pair_comparisons\":"<<comparisons<<",\"positive_diagonal_controls\":"<<positive_controls<<",\"obstruction_representations\":[";
 for(int k=0;k<4;k++){if(k)std::cout<<",";std::cout<<counts[k];}
 std::cout<<"],\"first_witnesses\":[";
 for(int k=0;k<4;k++){
   if(k)std::cout<<",";
   if(witnesses[k][0]<0){std::cout<<"null";continue;}
   std::cout<<"{\"u\":[";for(int j=0;j<n;j++){if(j)std::cout<<",";std::cout<<e.vectors[witnesses[k][0]][j];}
   std::cout<<"],\"v\":[";for(int j=0;j<n;j++){if(j)std::cout<<",";std::cout<<witnesses[k][2]*e.vectors[witnesses[k][1]][j];}
   std::cout<<"]}";
 }
 std::cout<<"],\"enumeration_and_preparation_seconds\":"<<std::chrono::duration<double>(enum_done-started).count()<<",\"binary_check_seconds\":"<<std::chrono::duration<double>(ended-enum_done).count()<<",\"total_seconds\":"<<std::chrono::duration<double>(ended-started).count()<<"}\n";
 return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
