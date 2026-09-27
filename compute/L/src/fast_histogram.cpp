// Exact integer Fincke--Pohst enumeration. GMP is used only for setup.
// At level i, rem = remaining * D[i+1], c = D[i+1] sum U[i,j]x[j],
// where D[k] is the leading k by k determinant. The next remainder is
// (rem*D[i]-(D[i+1]*x[i]+c)^2)/D[i+1]. Every division is checked exact.
// All 128-bit arithmetic is certified safe from coordinate and minor bounds
// before entering the recursion. Unsupported inputs fail, never fall back.
#include <gmpxx.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <cstdint>
static_assert(sizeof(long)>=8,"64-bit long required");
using Q=mpq_class;using Z=mpz_class;using I=__int128_t;using U=__uint128_t;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
unsigned long long isqrt(U n){
 unsigned long long k=(unsigned long long)std::sqrt((long double)n);
 while(U(k)*k>n)--k;
 while(U(k+1)*(k+1)<=n)++k;
 return k;
}
I floordiv(I a,I b){I q=a/b,r=a%b;return q-(r<0);}
I ceildiv(I a,I b){return -floordiv(-a,b);}
struct Enum{
 int n,bound;std::vector<long>D,cap,x;std::vector<std::vector<long>>C;
 std::vector<unsigned long long>hist;unsigned long long nodes=0,pairs=0;
 Enum(std::istream&in){
  in>>n>>bound;need(bool(in)&&n>0&&n<=100&&bound>0&&bound<=1000000,"header");
  std::vector<std::vector<Q>>a(n,std::vector<Q>(n)),u=a,inv=a;std::vector<Q>d(n);
  for(auto&r:a)for(auto&v:r){in>>v;v.canonicalize();need(v.get_den()==1,"integral input");}need(bool(in),"input complete");
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)need(a[i][j]==a[j][i],"symmetric");
  for(int i=0;i<n;i++){d[i]=a[i][i];for(int k=0;k<i;k++)d[i]-=d[k]*u[k][i]*u[k][i];need(d[i]>0,"positive definite");u[i][i]=1;for(int j=i+1;j<n;j++){Q v=a[i][j];for(int k=0;k<i;k++)v-=d[k]*u[k][i]*u[k][j];u[i][j]=v/d[i];}}
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){Q v=0;for(int k=0;k<n;k++)v+=d[k]*u[k][i]*u[k][j];need(v==a[i][j],"exact LDL identity");}
  auto aa=a;for(int i=0;i<n;i++)inv[i][i]=1;
  for(int i=0;i<n;i++){Q p=aa[i][i];need(p!=0,"inverse pivot");for(int j=0;j<n;j++){aa[i][j]/=p;inv[i][j]/=p;}for(int k=0;k<n;k++)if(k!=i){Q q=aa[k][i];for(int j=0;j<n;j++){aa[k][j]-=q*aa[i][j];inv[k][j]-=q*inv[i][j];}}}
  D={1};cap.resize(n);x.resize(n);C.assign(n,std::vector<long>(n));hist.resize(bound+1);Q det=1;Z safety=Z(1)<<62;
  for(int i=0;i<n;i++){
   det*=d[i];need(det.get_den()==1&&det>0&&det<safety,"minor range");D.push_back(det.get_num().get_si());
   Q lim=bound*inv[i][i];Z fl=lim.get_num()/lim.get_den(),rt;mpz_sqrt(rt.get_mpz_t(),fl.get_mpz_t());need(rt<safety,"coordinate bound");cap[i]=rt.get_si();
   for(int j=i+1;j<n;j++){Q v=det*u[i][j];need(v.get_den()==1&&abs(v)<safety,"scaled coefficient range");C[i][j]=v.get_num().get_si();}
  }
  for(int i=0;i<n;i++){
   Z cm=0;for(int j=i+1;j<n;j++)cm+=abs(Z(C[i][j]))*Z(cap[j]);
   need(cm+Z(D[i+1])*Z(cap[i])<safety,"linear expression bound");
   need(Z(bound)*Z(D[i])*Z(D[i+1])<(Z(1)<<124),"squared radius bound");
  }
 }
 void rec(int i,I rem,bool nz){
  if(i<0){if(nz){need(rem>=0&&rem<bound,"leaf norm");hist[bound-(int)rem]+=2;pairs++;}return;}
  long long c=0;for(int j=i+1;j<n;j++)c+=C[i][j]*x[j];
  I rad=rem*D[i];need(rad>=0,"nonnegative radius");long long k=isqrt(U(rad));
  I lo=ceildiv(-I(k)-c,D[i+1]),hi=floordiv(I(k)-c,D[i+1]);if(!nz&&lo<0)lo=0;
  if(lo>hi)return;need(lo>=-cap[i]&&hi<=cap[i],"coordinate bound respected");
  for(long long q=(long long)lo;q<=(long long)hi;q++){
   x[i]=q;nodes++;I v=I(D[i+1])*q+c,rr=rad-v*v;need(rr>=0&&rr%D[i+1]==0,"exact nonnegative next remainder");rec(i-1,rr/D[i+1],nz||q!=0);
  }
 }
};
int main(int argc,char**argv){try{
 need(argc==2||argc==3,"usage fast_histogram [--batch] file");bool batch=argc==3;need(!batch||std::string(argv[1])=="--batch","batch flag");std::ifstream in(argv[batch?2:1]);need(bool(in),"open input");int idx=0;
 while(in>>std::ws && in.peek()!=EOF){auto t=std::chrono::steady_clock::now();Enum e(in);e.rec(e.n-1,I(e.bound)*e.D[e.n],false);std::cout<<"{\"form\":"<<++idx<<",\"dimension\":"<<e.n<<",\"bound\":"<<e.bound<<",\"complete\":true,\"pairs\":"<<e.pairs<<",\"nodes\":"<<e.nodes<<",\"histogram\":[";for(size_t j=0;j<e.hist.size();j++){if(j)std::cout<<',';std::cout<<e.hist[j];}std::cout<<"],\"wall_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<"}"<<std::endl;if(!batch)break;}
 }catch(std::exception const&e){std::cerr<<e.what()<<'\n';return 1;}}
