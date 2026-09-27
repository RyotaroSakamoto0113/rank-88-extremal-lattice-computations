#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <numeric>
#include <chrono>
#include <stdexcept>
using Vec=std::array<int,22>; using Mat=std::array<Vec,22>;
struct Hash {size_t operator()(Vec const& a)const{size_t s=1469598103934665603ULL;for(int x:a){s^=size_t(x);s*=1099511628211ULL;}return s;}};
Vec canon(Vec v){for(int x:v)if(x){if(x<0)for(int &y:v)y=-y;break;}return v;}
Vec mul(Mat const&a,Vec const&v){Vec r{};for(int i=0;i<22;++i)for(int j=0;j<22;++j)r[i]+=a[i][j]*v[j];return r;}
int dot(Vec const&a,Vec const&b){int r=0;for(int i=0;i<22;++i)r+=a[i]*b[i];return r;}
int main(){auto t=std::chrono::steady_clock::now();std::ifstream in("output/rank4_audit/M_integer_data.txt");int n,ng;in>>n;if(n!=22)return 1;Mat G,W;for(auto&m:{&G,&W})for(auto&r:*m)for(int&x:r)in>>x;in>>ng;std::vector<Mat> gen(ng);for(auto&m:gen)for(auto&r:m)for(int&x:r)in>>x;if(!in)return 2;
 std::ifstream sf("output/rank4_audit/M_shell12.txt");int nv;sf>>nv>>n;std::vector<Vec> vs(nv);std::vector<int> norms(nv),par(nv);std::iota(par.begin(),par.end(),0);auto find=[&](int x){while(x!=par[x])x=par[x]=par[par[x]];return x;};std::unordered_map<Vec,int,Hash> ids;ids.reserve(nv*2);std::map<int,int> shell,orbits;for(int i=0;i<nv;++i){for(int&x:vs[i])sf>>x;vs[i]=canon(vs[i]);auto [it,ok]=ids.emplace(vs[i],i);if(!ok)throw std::runtime_error("duplicate");int q=dot(vs[i],mul(G,vs[i]));if(q<=0||q>24||q%2)throw std::runtime_error("norm");norms[i]=q/2;++shell[q/2];}if(!sf)return 3;for(int i=0;i<nv;++i)for(auto&m:gen){Vec u=canon(mul(m,vs[i]));auto it=ids.find(u);if(it==ids.end())throw std::runtime_error("missing image");int a=find(i),b=find(it->second);if(a!=b)par[b]=a;}
 std::ofstream reps("output/rank4_audit/M_shell12_orbit_reps.txt");for(int i=0;i<nv;++i)if(find(i)==i){++orbits[norms[i]];reps<<norms[i];for(int x:vs[i])reps<<' '<<x;reps<<'\n';}
 std::cout<<"{\"pair_count\":"<<nv<<",\"signed_count\":"<<2*nv<<",\"unique_valid\":true,\"generator_closed\":true,\"shells\":[";bool sep=0;for(auto [q,c]:shell){if(sep)std::cout<<',';sep=1;std::cout<<"{\"q\":"<<q<<",\"pairs\":"<<c<<",\"orbits\":"<<orbits[q]<<"}";}std::cout<<"],\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<"}\n";
}
