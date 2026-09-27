#include <gmpxx.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using Q=mpq_class; using Z=mpz_class;
using Clock=std::chrono::steady_clock;
struct Affine {
 int n; Q radius; std::vector<std::vector<Q>> a,u; std::vector<Q>d,c;
 std::vector<Z>x; std::vector<std::vector<Z>> solutions;
 unsigned long long nodes=0; double limit; Clock::time_point start;
 Affine(std::istream&in,double lim):limit(lim){
  in>>n>>radius; radius.canonicalize();
  if(!in||n<1||n>100)throw std::runtime_error("invalid affine header");
  a.assign(n,std::vector<Q>(n));u=a;d.resize(n);c.resize(n);x.resize(n);
  for(auto&v:c){in>>v;v.canonicalize();}
  for(auto&row:a)for(auto&v:row){in>>v;v.canonicalize();}
  if(!in)throw std::runtime_error("truncated affine input");
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[i][j]!=a[j][i])throw std::runtime_error("asymmetric matrix");
  for(int i=0;i<n;i++){
   d[i]=a[i][i];for(int k=0;k<i;k++)d[i]-=d[k]*u[k][i]*u[k][i];
   if(d[i]<=0)throw std::runtime_error("nonpositive pivot");u[i][i]=1;
   for(int j=i+1;j<n;j++){Q t=a[i][j];for(int k=0;k<i;k++)t-=d[k]*u[k][i]*u[k][j];u[i][j]=t/d[i];}
  }
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){Q t=0;for(int k=0;k<n;k++)t+=d[k]*u[k][i]*u[k][j];if(t!=a[i][j])throw std::runtime_error("LDL reconstruction failed");}
 }
 void rec(int i,const Q&left){
  if((nodes&8191)==0 && std::chrono::duration<double>(Clock::now()-start).count()>limit)throw std::runtime_error("affine time limit: incomplete");
  if(i<0){
   Q check=0;for(int j=0;j<n;j++)for(int k=0;k<n;k++)check+=(Q(x[j])-c[j])*a[j][k]*(Q(x[k])-c[k]);
   if(check>radius)throw std::runtime_error("invalid affine witness");solutions.push_back(x);return;
  }
  Q offset=-c[i];for(int j=i+1;j<n;j++)offset+=u[i][j]*(Q(x[j])-c[j]);
  Q rad=left/d[i]; Z den=offset.get_den(),num=offset.get_num();
  Z scale=rad.get_num()*den*den, div=rad.get_den(),floor;
  mpz_fdiv_q(floor.get_mpz_t(),scale.get_mpz_t(),div.get_mpz_t());
  Z rt;mpz_sqrt(rt.get_mpz_t(),floor.get_mpz_t());
  Z lo,hi,v=-rt-num;mpz_cdiv_q(lo.get_mpz_t(),v.get_mpz_t(),den.get_mpz_t());
  v=rt-num;mpz_fdiv_q(hi.get_mpz_t(),v.get_mpz_t(),den.get_mpz_t());
  for(Z t=lo;t<=hi;t++){
   x[i]=t;nodes++;Q y=Q(t)+offset,rem=left-d[i]*y*y;
   if(rem<0)throw std::runtime_error("affine interval error");rec(i-1,rem);
  }
 }
 void run(){start=Clock::now();if(radius>=0)rec(n-1,radius);}
};
int main(int argc,char**argv){try{
 if(argc<2)throw std::runtime_error("usage: exact_affine_enum input [limit_seconds]");
 std::ifstream in(argv[1]);if(!in)throw std::runtime_error("cannot open input");
 double limit=argc>2?std::stod(argv[2]):60;int id=0;
 while(in>>std::ws && in.peek()!=EOF){
  auto start=Clock::now();Affine e(in,limit);e.run();
  std::cout<<"{\"id\":"<<++id<<",\"complete\":true,\"dimension\":"<<e.n<<",\"count\":"<<e.solutions.size()<<",\"nodes\":"<<e.nodes<<",\"seconds\":"<<std::chrono::duration<double>(Clock::now()-start).count()<<",\"solutions\":[";
  for(size_t j=0;j<e.solutions.size();j++){if(j)std::cout<<",";std::cout<<"[";for(int k=0;k<e.n;k++){if(k)std::cout<<",";std::cout<<e.solutions[j][k];}std::cout<<"]";}
  std::cout<<"]}\n";
 }
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
