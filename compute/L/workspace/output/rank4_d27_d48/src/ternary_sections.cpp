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
struct Enum8{std::array<std::array<Q,8>,8> u{};std::array<Q,8>d{};V8 x{};std::vector<SV> vs;const M8&G;int bound=12;unsigned long nodes=0;
 Enum8(const M8&a,int bb=12):G(a),bound(bb){for(int i=0;i<8;++i){d[i]=a[i][i];for(int k=0;k<i;++k)d[i]-=d[k]*u[k][i]*u[k][i];check(d[i]>0,"not positive definite");u[i][i]=1;for(int j=i+1;j<8;++j){Q v=a[i][j];for(int k=0;k<i;++k)v-=d[k]*u[k][i]*u[k][j];u[i][j]=v/d[i];}}for(int i=0;i<8;++i)for(int j=0;j<8;++j){Q v=0;for(int k=0;k<8;++k)v+=d[k]*u[k][i]*u[k][j];check(v==a[i][j],"LDL identity");}}
 void rec(int i,Q remaining,bool nz){if(i<0){if(nz){SV v;v.x=x;for(int j=0;j<8;++j)for(int k=0;k<8;++k)v.g[j]+=x[k]*G[k][j];for(int j=0;j<4;++j){v.gw[j]=v.g[j+4];v.gw[j+4]=-6*v.g[j]+v.g[j+4];}long q=0;for(int j=0;j<8;++j)q+=x[j]*v.g[j];check(q>0&&q<=2*bound&&q%2==0,"source vector norm");v.q=q/2;vs.push_back(v);}return;}Q c=0;for(int j=i+1;j<8;++j)c+=u[i][j]*x[j];Q r=remaining/d[i];Z aa=c.get_num(),bb=c.get_den(),sc=r.get_num()*bb*bb,den=r.get_den(),ir,k;mpz_fdiv_q(ir.get_mpz_t(),sc.get_mpz_t(),den.get_mpz_t());mpz_sqrt(k.get_mpz_t(),ir.get_mpz_t());Z lo,hi,nn=-k-aa;mpz_cdiv_q(lo.get_mpz_t(),nn.get_mpz_t(),bb.get_mpz_t());nn=k-aa;mpz_fdiv_q(hi.get_mpz_t(),nn.get_mpz_t(),bb.get_mpz_t());if(!nz&&lo<0)lo=0;check(lo.fits_slong_p()&&hi.fits_slong_p(),"coordinate overflow");long last=hi.get_si();for(long t=lo.get_si();t<=last;++t){check(std::abs(t)<1000000,"unexpected large source coordinate");x[i]=t;++nodes;Q v=c+t,rem=remaining-d[i]*v*v;check(rem>=0,"radius");rec(i-1,rem,nz||t!=0);if(t==last)break;}}
};
Pair h(SV const&u,SV const&v){long s=0,t=0;for(int i=0;i<8;++i){s+=u.g[i]*v.x[i];t+=u.gw[i]*v.x[i];}check((12*s-t)%23==0&&(2*t-s)%23==0,"nonintegral source h");return{(12*s-t)/23,(2*t-s)/23};}
M8 gram(H4 const&H){M8 g{};for(int i=0;i<4;++i)for(int j=0;j<4;++j){auto z=H[i][j];auto y=H[j][i];check(y.a==z.a+z.b&&y.b==-z.b,"not Hermitian");check(std::abs(z.a)<1000000000&&std::abs(z.b)<1000000000,"large Gram coefficient");g[i][j]=2*z.a+z.b;g[i][j+4]=z.a+12*z.b;g[i+4][j]=z.a-11*z.b;g[i+4][j+4]=6*(2*z.a+z.b);}for(int i=0;i<8;++i)for(int j=0;j<8;++j)check(g[i][j]==g[j][i],"trace symmetry");return g;}
int rank(std::vector<V8> const&v){int n=v.size();std::vector<std::array<Q,8>>a(n);for(int i=0;i<n;++i)for(int j=0;j<8;++j)a[i][j]=v[i][j];int k=0;for(int j=0;j<8&&k<n;++j){int p=k;while(p<n&&a[p][j]==0)++p;if(p==n)continue;std::swap(a[k],a[p]);Q z=a[k][j];for(int c=j;c<8;++c)a[k][c]/=z;for(int r=k+1;r<n;++r){z=a[r][j];for(int c=j;c<8;++c)a[r][c]-=z*a[k][c];}++k;}return k;}
V8 wmul(V8 const&x){V8 r;for(int i=0;i<4;++i){r[i]=-6*x[i+4];r[i+4]=x[i]+x[i+4];}return r;}
void vecout(V8 const&v){std::cout<<'[';for(int i=0;i<8;++i){if(i)std::cout<<',';std::cout<<v[i];}std::cout<<']';}

Pair mul(Pair x,Pair y){return{x.a*y.a-6*x.b*y.b,x.a*y.b+x.b*y.a+x.b*y.b};}
Pair conj(Pair z){return{z.a+z.b,-z.b};}
long norm(Pair z){return z.a*z.a+z.a*z.b+6*z.b*z.b;}
long det3(int a,int b,int c,Pair x,Pair y,Pair z){Pair p=mul(mul(x,z),conj(y));return long(a)*b*c-a*norm(z)-b*norm(y)-c*norm(x)+2*p.a+p.b;}
void hmatrixout(std::vector<int>const& ix,std::vector<SV>const&v){std::cout<<'[';for(size_t i=0;i<ix.size();++i){if(i)std::cout<<',';std::cout<<'[';for(size_t j=0;j<ix.size();++j){if(j)std::cout<<',';auto z=h(v[ix[i]],v[ix[j]]);std::cout<<'['<<z.a<<','<<z.b<<']';}std::cout<<']';}std::cout<<']';}
int main(int argc,char**argv){try{
 check(argc>=2&&argc<=4,"usage ternary_sections flat_input [bound=12] [mode=0 cached / 1 direct]");
 int bound=argc>2?std::stoi(argv[2]):12, mode=argc>3?std::stoi(argv[3]):0;
 check(bound>=1&&bound<=30,"bound range");check(mode==0||mode==1,"mode range");
 auto start=Clock::now();std::ifstream in(argv[1]);check(bool(in),"input missing");long id;int total=0;std::map<std::string,int>num;unsigned long pairtests=0,tripletests=0,nodes=0;double enumsec=0,pairsec=0,triplesec=0;
 while(in>>id){auto t0=Clock::now();H4 H;for(auto&r:H)for(auto&z:r)in>>z.a>>z.b;check(bool(in),"incomplete source");M8 G=gram(H);Enum8 e(G,bound);e.rec(7,Q(2*bound),false);auto&v=e.vs;
 std::sort(v.begin(),v.end(),[](auto&a,auto&b){if(a.q!=b.q)return a.q<b.q;return a.x<b.x;});auto t1=Clock::now();
 std::string status="survives";std::vector<int>witness;long det=0;unsigned long np=0,nt=0;
 for(int i=0;i<(int)v.size();++i)if(v[i].q<6){status="det1_lt6";witness={i};det=v[i].q;break;}
 int n=v.size();std::vector<Pair> table(n*n);std::vector<long> d2(n*n);
 if(status=="survives")for(int i=0;i<n;++i){for(int j=i+1;j<n;++j){Pair z=h(v[i],v[j]);table[i*n+j]=z;long d=long(v[i].q)*v[j].q-norm(z);d2[i*n+j]=d;++np;check(d>=0,"negative binary det");if(d>0&&d<16){status="det2_lt16";witness={i,j};det=d;break;}}if(!witness.empty())break;}
 auto t2=Clock::now();
 if(status=="survives")for(int i=0;i<n;++i){for(int j=i+1;j<n;++j){if(!d2[i*n+j])continue;for(int k=j+1;k<n;++k){if(!d2[i*n+k]||!d2[j*n+k])continue;Pair a,b,c;if(mode==0){a=table[i*n+j];b=table[i*n+k];c=table[j*n+k];}else{a=h(v[i],v[j]);b=h(v[i],v[k]);c=h(v[j],v[k]);}long d=det3(v[i].q,v[j].q,v[k].q,a,b,c);++nt;check(d>=0,"negative ternary det");if(d>0&&d<27){status="det3_lt27";witness={i,j,k};det=d;break;}}if(!witness.empty())break;}if(!witness.empty())break;}
 auto t3=Clock::now();double se=std::chrono::duration<double>(t1-t0).count(),sp=std::chrono::duration<double>(t2-t1).count(),st=std::chrono::duration<double>(t3-t2).count();
 std::cout<<"{\"id\":"<<id<<",\"status\":\""<<status<<"\",\"short_pairs\":"<<v.size()<<",\"enum_nodes\":"<<e.nodes<<",\"pair_tests\":"<<np<<",\"triple_tests\":"<<nt<<",\"enum_seconds\":"<<se<<",\"pair_seconds\":"<<sp<<",\"triple_seconds\":"<<st;
 if(!witness.empty()){std::cout<<",\"witness\":[";for(size_t k=0;k<witness.size();++k){if(k)std::cout<<',';vecout(v[witness[k]].x);}std::cout<<"],\"witness_Gram\":";hmatrixout(witness,v);std::cout<<",\"det\":"<<det;}
 std::cout<<"}\n";++total;++num[status];nodes+=e.nodes;pairtests+=np;tripletests+=nt;enumsec+=se;pairsec+=sp;triplesec+=st;
 }
 double wall=std::chrono::duration<double>(Clock::now()-start).count();std::cerr<<"{\"total\":"<<total<<",\"bound\":"<<bound<<",\"mode\":"<<mode<<",\"enum_nodes\":"<<nodes<<",\"pair_tests\":"<<pairtests<<",\"triple_tests\":"<<tripletests<<",\"enum_seconds\":"<<enumsec<<",\"pair_seconds\":"<<pairsec<<",\"triple_seconds\":"<<triplesec<<",\"total_seconds\":"<<wall<<",\"counts\":{";bool first=true;for(auto&[s,n]:num){if(!first)std::cerr<<',';first=false;std::cerr<<'\"'<<s<<"\":"<<n;}std::cerr<<"}}\n";return 0;
 }catch(std::exception const&e){std::cerr<<"ERROR "<<e.what()<<'\n';return 1;}}
