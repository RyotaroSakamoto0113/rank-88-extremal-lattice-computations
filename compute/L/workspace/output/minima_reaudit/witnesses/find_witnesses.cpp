#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <stdexcept>
#include <string>
using V=std::array<long,22>;
using M=std::array<V,22>;
struct P {long a,b;};
bool eq(P x,P y){return x.a==y.a&&x.b==y.b;}
void need(bool ok,const char*s){if(!ok)throw std::runtime_error(s);}
struct S {V x{},g{},gw{};int q;};
M G{},W{};std::vector<S> v;
P ip(S const&x,S const&y){long s=0,t=0;for(int i=0;i<22;++i){s+=x.g[i]*y.x[i];t+=x.gw[i]*y.x[i];}need((12*s-t)%23==0&&(2*t-s)%23==0,"nonintegral pairing");return {(12*s-t)/23,(2*t-s)/23};}
P ip(int i,int j){P z=ip(v[std::abs(i)-1],v[std::abs(j)-1]);if((i<0)!=(j<0))return{-z.a,-z.b};return z;}
std::vector<int> list(int a,int q,P z){std::vector<int> out;for(int i=1;i<=(int)v.size();++i)if(v[i-1].q==q){P h=ip(a,i);if(eq(h,z))out.push_back(i);if(eq({-h.a,-h.b},z))out.push_back(-i);}return out;}
void output(int r,std::vector<int>const&ids){std::cout<<"{\"rank\":"<<r<<",\"shell_ids\":[";for(int j=0;j<r;++j){if(j)std::cout<<',';std::cout<<ids[j];}std::cout<<"],\"coordinates\":[";for(int j=0;j<r;++j){if(j)std::cout<<',';std::cout<<'[';for(int i=0;i<22;++i){if(i)std::cout<<',';std::cout<<(ids[j]<0?-1:1)*v[std::abs(ids[j])-1].x[i];}std::cout<<']';}std::cout<<"],\"gram\":[";for(int j=0;j<r;++j){if(j)std::cout<<',';std::cout<<'[';for(int k=0;k<r;++k){if(k)std::cout<<',';auto h=ip(ids[j],ids[k]);std::cout<<'['<<h.a<<','<<h.b<<']';}std::cout<<']';}std::cout<<"]}\n";std::cout.flush();}
int main(int argc,char**argv){try{need(argc==4,"usage find_witnesses M_integer_data.txt M_shell12.txt M_shell12_orbit_reps.txt");auto t=std::chrono::steady_clock::now();std::ifstream f(argv[1]);int n;f>>n;need(n==22,"ambient dimension");for(auto m:{&G,&W})for(auto&r:*m)for(auto&x:r)f>>x;need(bool(f),"ambient read");std::ifstream s(argv[2]);int count;s>>count>>n;need(n==22,"shell dimension");v.resize(count);int anchor=0;for(int k=0;k<count;++k){auto&a=v[k];for(long&x:a.x)s>>x;for(int j=0;j<22;++j)for(int i=0;i<22;++i)a.g[j]+=a.x[i]*G[i][j];for(int j=0;j<22;++j)for(int i=0;i<22;++i)a.gw[j]+=a.g[i]*W[i][j];auto h=ip(a,a);need(h.b==0&&h.a>0&&h.a<=12,"shell norm");a.q=h.a;if(!anchor&&a.q==6)anchor=k+1;}need(bool(s)&&anchor,"shell read");output(1,{anchor});std::cerr<<"input_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<'\n';
bool found=false;std::ifstream reps(argv[3]);int rq;long binary_tests=0;while(reps>>rq&&!found){V u;for(long&x:u)reps>>x;int a=0;for(int i=1;i<=count&&!a;++i)if(v[i-1].x==u)a=i;need(a,"orbit representative absent from shell");for(int i=1;i<=count&&!found;++i){++binary_tests;auto h=ip(a,i);if(rq*v[i-1].q-h.a*h.a-h.a*h.b-6*h.b*h.b==16){output(2,{a,i});found=true;}}}std::cerr<<"d2_pair_tests "<<binary_tests<<'\n';need(found,"no d2 witness found");
auto b3=list(anchor,9,{-3,0}),c3=list(anchor,12,{0,-3});long triples=0;found=false;for(int b:b3){for(int c:c3){++triples;if(eq(ip(b,c),{3,0})){output(3,{anchor,b,c});found=true;break;}}if(found)break;}std::cerr<<"d3_first_extension "<<b3.size()<<" d3_second_extension "<<c3.size()<<" d3_pair_tests "<<triples<<'\n';need(found,"no d3 witness found");
auto b4=list(anchor,6,{-2,0}),c4=list(anchor,8,{-2,2}),d4=list(anchor,10,{-2,-2});long pairs=0,tests=0;found=false;for(int b:b4){std::vector<int> last;for(int d:d4)if(eq(ip(b,d),{2,0}))last.push_back(d);for(int c:c4){++pairs;if(!eq(ip(b,c),{0,0}))continue;for(int d:last){++tests;if(eq(ip(c,d),{-2,0})){output(4,{anchor,b,c,d});found=true;break;}}if(found)break;}if(found)break;}std::cerr<<"d4_extensions "<<b4.size()<<' '<<c4.size()<<' '<<d4.size()<<" d4_pair_tests "<<pairs<<" d4_last_tests "<<tests<<'\n';need(found,"no d4 witness found");std::cerr<<"total_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<'\n';return 0;}catch(std::exception const&e){std::cerr<<"ERROR "<<e.what()<<'\n';return 1;}}
