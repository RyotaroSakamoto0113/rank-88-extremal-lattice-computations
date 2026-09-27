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
 std::vector<long>x,first; unsigned long long nodes=0,pairs=0; bool complete=true; std::vector<unsigned long long> histogram;
 double seconds_limit; Clock::time_point started;
 Enumerator(std::istream&in,double limit):seconds_limit(limit){
  in>>n>>bound; if(!in||n<1||n>100||bound<0)throw std::runtime_error("invalid header");bound.canonicalize(); if(bound.get_den()!=1 || bound>1000000)throw std::runtime_error("histogram bound"); histogram.resize(bound.get_num().get_ui()+1);
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
  if(i<0){if(nonzero){Q nn=bound-remaining;if(nn.get_den()!=1)throw std::runtime_error("noninteger norm"); histogram.at(nn.get_num().get_ui())++;pairs++;if(first.empty())first=x;}return;}
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
bool process(std::istream&in,const std::string&label,double limit){
  auto begin=Clock::now();
  Enumerator e(in,limit);auto reduced=Clock::now();e.started=reduced;
  try{e.rec(e.n-1,e.bound,false);}catch(const std::runtime_error&ex){if(e.complete)throw;}
  auto end=Clock::now();
  if(!e.first.empty()) { Q norm=0;for(int i=0;i<e.n;i++)for(int j=0;j<e.n;j++)norm+=e.first[i]*e.a[i][j]*e.first[j];if(norm<=0||norm>e.bound)throw std::runtime_error("witness check failed"); }
  std::cout<<"{\"file\":\""<<label<<"\",\"dimension\":"<<e.n<<",\"bound\":\""<<e.bound<<"\",\"complete\":"<<(e.complete?"true":"false")<<",\"pairs\":"<<e.pairs<<",\"signed_vectors\":"<<2*e.pairs<<",\"nodes\":"<<e.nodes<<",\"ldl_seconds\":"<<std::chrono::duration<double>(reduced-begin).count()<<",\"enumeration_seconds\":"<<std::chrono::duration<double>(end-reduced).count()<<",\"first\":[";
  for(size_t i=0;i<e.first.size();i++){if(i)std::cout<<",";std::cout<<e.first[i];}std::cout<<"],\"histogram\":[";for(size_t j=0;j<e.histogram.size();j++){if(j)std::cout<<",";std::cout<<2*e.histogram[j];}std::cout<<"]}\n";
  return e.complete;
}
int main(int argc,char**argv){
 try{
  if(argc<2)throw std::runtime_error("usage: exact_enum [--batch] input [seconds_limit]");
  bool batch=std::string(argv[1])=="--batch";int pos=batch?2:1;
  if(argc<=pos)throw std::runtime_error("missing file");
  std::ifstream in(argv[pos]);if(!in)throw std::runtime_error("cannot open input");
  double limit=argc>pos+1?std::stod(argv[pos+1]):60;
  if(!batch)return process(in,argv[pos],limit)?0:2;
  unsigned long long index=0;
  while(in>>std::ws && in.peek()!=EOF){index++;if(!process(in,std::to_string(index),limit))return 2;}
  return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
}
