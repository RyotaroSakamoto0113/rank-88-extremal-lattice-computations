#include <gmpxx.h>
#include <array>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>
using Q=mpq_class;using Z=mpz_class;using Clock=std::chrono::steady_clock;using V8=std::array<long,8>;using V22=std::array<long,22>;using M8=std::array<V8,8>;using M22=std::array<V22,22>;
struct Pair{long a=0,b=0;};using H4=std::array<std::array<Pair,4>,4>;
struct SV{V8 x{},g{},gw{};int q;};
void check(bool b,const char*s){if(!b)throw std::runtime_error(s);}
struct Enum8{std::array<std::array<Q,8>,8> u{};std::array<Q,8>d{};V8 x{};std::vector<SV> vs;const M8&G;unsigned long nodes=0;
 Enum8(const M8&a):G(a){for(int i=0;i<8;++i){d[i]=a[i][i];for(int k=0;k<i;++k)d[i]-=d[k]*u[k][i]*u[k][i];check(d[i]>0,"not positive definite");u[i][i]=1;for(int j=i+1;j<8;++j){Q v=a[i][j];for(int k=0;k<i;++k)v-=d[k]*u[k][i]*u[k][j];u[i][j]=v/d[i];}}for(int i=0;i<8;++i)for(int j=0;j<8;++j){Q v=0;for(int k=0;k<8;++k)v+=d[k]*u[k][i]*u[k][j];check(v==a[i][j],"LDL identity");}}
 void rec(int i,Q remaining,bool nz){if(i<0){if(nz){SV v;v.x=x;for(int j=0;j<8;++j)for(int k=0;k<8;++k)v.g[j]+=x[k]*G[k][j];for(int j=0;j<4;++j){v.gw[j]=v.g[j+4];v.gw[j+4]=-6*v.g[j]+v.g[j+4];}long q=0;for(int j=0;j<8;++j)q+=x[j]*v.g[j];check(q>0&&q<=24&&q%2==0,"source vector norm");v.q=q/2;vs.push_back(v);}return;}Q c=0;for(int j=i+1;j<8;++j)c+=u[i][j]*x[j];Q r=remaining/d[i];Z aa=c.get_num(),bb=c.get_den(),sc=r.get_num()*bb*bb,den=r.get_den(),ir,k;mpz_fdiv_q(ir.get_mpz_t(),sc.get_mpz_t(),den.get_mpz_t());mpz_sqrt(k.get_mpz_t(),ir.get_mpz_t());Z lo,hi,nn=-k-aa;mpz_cdiv_q(lo.get_mpz_t(),nn.get_mpz_t(),bb.get_mpz_t());nn=k-aa;mpz_fdiv_q(hi.get_mpz_t(),nn.get_mpz_t(),bb.get_mpz_t());if(!nz&&lo<0)lo=0;check(lo.fits_slong_p()&&hi.fits_slong_p(),"coordinate overflow");long last=hi.get_si();for(long t=lo.get_si();t<=last;++t){check(std::abs(t)<1000000,"unexpected large source coordinate");x[i]=t;++nodes;Q v=c+t,rem=remaining-d[i]*v*v;check(rem>=0,"radius");rec(i-1,rem,nz||t!=0);if(t==last)break;}}
};
Pair h(SV const&u,SV const&v){long s=0,t=0;for(int i=0;i<8;++i){s+=u.g[i]*v.x[i];t+=u.gw[i]*v.x[i];}check((12*s-t)%23==0&&(2*t-s)%23==0,"nonintegral source h");return{(12*s-t)/23,(2*t-s)/23};}
constexpr int AS=65,BS=33;constexpr long SZ=13*13*AS*BS;
long key(int q,int r,Pair z){if(q>r){std::swap(q,r);z={z.a+z.b,-z.b};}check(q>=0&&r<=12&&std::abs(z.a)<=32&&std::abs(z.b)<=16,"pair key range");return((q*13+r)*AS+z.a+32)*BS+z.b+16;}
struct Ambient{M22 G{},W{};std::vector<V22>vs;std::vector<int>qn;std::vector<unsigned> counts;std::array<int,13>shell{};Ambient():counts(SZ){std::ifstream f("output/rank4_audit/M_integer_data.txt");int n;f>>n;check(n==22,"ambient dimension");for(auto m:{&G,&W})for(auto&r:*m)for(auto&x:r)f>>x;check(bool(f),"ambient read");std::ifstream s("output/rank4_audit/M_shell12.txt");int nv;s>>nv>>n;check(n==22,"shell dimension");vs.resize(nv);qn.resize(nv);for(int k=0;k<nv;++k){for(long&x:vs[k])s>>x;long norm=0;for(int i=0;i<22;++i)for(int j=0;j<22;++j)norm+=vs[k][i]*G[i][j]*vs[k][j];check(norm>0&&norm<=24&&norm%2==0,"ambient norm");qn[k]=norm/2;++shell[norm/2];}check(bool(s),"shell read");}
 void table(){std::ifstream f("output/rank4_audit/M_shell12_orbit_reps.txt");int q;unsigned long tests=0;while(f>>q){V22 u,ug{},ugw{};for(long&x:u)f>>x;for(int j=0;j<22;++j)for(int i=0;i<22;++i)ug[j]+=u[i]*G[i][j];for(int j=0;j<22;++j)for(int i=0;i<22;++i)ugw[j]+=ug[i]*W[i][j];for(size_t k=0;k<vs.size();++k)if(qn[k]>=q){long s=0,t=0;for(int j=0;j<22;++j){s+=ug[j]*vs[k][j];t+=ugw[j]*vs[k][j];}check((12*s-t)%23==0&&(2*t-s)%23==0,"nonintegral ambient h");Pair z{(12*s-t)/23,(2*t-s)/23};++counts[key(q,qn[k],z)];++counts[key(q,qn[k],{-z.a,-z.b})];tests+=2;}}std::cerr<<"ambient direct signed comparisons "<<tests<<"\n";}
};
M8 gram(H4 const&H){M8 g{};for(int i=0;i<4;++i)for(int j=0;j<4;++j){auto z=H[i][j];auto y=H[j][i];check(y.a==z.a+z.b&&y.b==-z.b,"not Hermitian");check(std::abs(z.a)<1000000000&&std::abs(z.b)<1000000000,"large Gram coefficient");g[i][j]=2*z.a+z.b;g[i][j+4]=z.a+12*z.b;g[i+4][j]=z.a-11*z.b;g[i+4][j+4]=6*(2*z.a+z.b);}for(int i=0;i<8;++i)for(int j=0;j<8;++j)check(g[i][j]==g[j][i],"trace symmetry");return g;}
int rank(std::vector<V8> const&v){int n=v.size();std::vector<std::array<Q,8>>a(n);for(int i=0;i<n;++i)for(int j=0;j<8;++j)a[i][j]=v[i][j];int k=0;for(int j=0;j<8&&k<n;++j){int p=k;while(p<n&&a[p][j]==0)++p;if(p==n)continue;std::swap(a[k],a[p]);Q z=a[k][j];for(int c=j;c<8;++c)a[k][c]/=z;for(int r=k+1;r<n;++r){z=a[r][j];for(int c=j;c<8;++c)a[r][c]-=z*a[k][c];}++k;}return k;}
V8 wmul(V8 const&x){V8 r;for(int i=0;i<4;++i){r[i]=-6*x[i+4];r[i+4]=x[i]+x[i+4];}return r;}
void vecout(V8 const&v){std::cout<<'[';for(int i=0;i<8;++i){if(i)std::cout<<',';std::cout<<v[i];}std::cout<<']';}
int main(int argc,char**argv){try{check(argc==2,"usage independent_embedding_sieve flat_input");auto t0=Clock::now();Ambient M;M.table();auto t1=Clock::now();std::ifstream in(argv[1]);long id;int total=0;std::map<std::string,int>num;while(in>>id){H4 H;for(auto&r:H)for(auto&z:r)in>>z.a>>z.b;check(bool(in),"incomplete source");M8 G=gram(H);Enum8 e(G);e.rec(7,Q(24),false);auto&v=e.vs;std::sort(v.begin(),v.end(),[](auto&a,auto&b){if(a.q!=b.q)return a.q<b.q;return a.x<b.x;});std::string status="survives";int ia=-1,ib=-1;Pair hh;long det=0;for(int i=0;i<(int)v.size();++i)if(!M.shell[v[i].q]){status="unary_absent";ia=i;break;}
 if(status=="survives")for(int i=0;i<(int)v.size();++i){for(int j=i+1;j<(int)v.size();++j){Pair z=h(v[i],v[j]);long d=long(v[i].q)*v[j].q-z.a*z.a-z.a*z.b-6*z.b*z.b;check(d>=0,"binary det negative");if(d>0&&d<16){status="det_lt16";ia=i;ib=j;hh=z;det=d;break;}}if(ia>=0)break;}
 if(status=="survives")for(int i=0;i<(int)v.size();++i){for(int j=i+1;j<(int)v.size();++j){Pair z=h(v[i],v[j]);long d=long(v[i].q)*v[j].q-z.a*z.a-z.a*z.b-6*z.b*z.b;if(d>0&&!M.counts[key(v[i].q,v[j].q,z)]){status="binary_absent";ia=i;ib=j;hh=z;det=d;break;}}if(ia>=0)break;}
 std::cout<<"{\"id\":"<<id<<",\"status\":\""<<status<<"\",\"short_pairs\":"<<v.size()<<",\"enum_nodes\":"<<e.nodes;
 if(ia>=0){std::cout<<",\"u\":";vecout(v[ia].x);std::cout<<",\"q_u\":"<<v[ia].q;}if(ib>=0){std::cout<<",\"v\":";vecout(v[ib].x);std::cout<<",\"q_v\":"<<v[ib].q<<",\"h\":["<<hh.a<<','<<hh.b<<"],\"det\":"<<det;}
 if(status=="survives"){std::vector<V8>span;std::vector<int>ix;for(int i=0;i<(int)v.size()&&ix.size()<4;++i){auto test=span;test.push_back(v[i].x);test.push_back(wmul(v[i].x));if(rank(test)>(int)span.size()){span=test;ix.push_back(i);}}std::cout<<",\"short_rank\":"<<ix.size()<<",\"frame\":[";for(size_t j=0;j<ix.size();++j){if(j)std::cout<<',';vecout(v[ix[j]].x);}std::cout<<']';}
 std::cout<<"}\n";++num[status];++total;}
 auto t2=Clock::now();std::cerr<<"total "<<total<<" ambient_sec "<<std::chrono::duration<double>(t1-t0).count()<<" candidate_sec "<<std::chrono::duration<double>(t2-t1).count()<<"\n";for(auto&[s,n]:num)std::cerr<<s<<' '<<n<<'\n';return 0;}catch(std::exception const&e){std::cerr<<"ERROR "<<e.what()<<'\n';return 1;}}
