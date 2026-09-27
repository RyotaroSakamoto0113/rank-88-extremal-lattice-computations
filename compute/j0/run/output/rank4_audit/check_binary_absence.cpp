// Independent direct inner-product check for the explicit forbidden types.
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using I=long long;using V=std::array<I,22>;using Mat=std::array<V,22>;
V mv(const Mat&A,const V&v){V o{};for(int i=0;i<22;i++)for(int j=0;j<22;j++)o[i]+=A[i][j]*v[j];return o;}
I dot(const V&u,const V&v){I r=0;for(int i=0;i<22;i++)r+=u[i]*v[i];return r;}
int main(){try{
 auto start=std::chrono::steady_clock::now();
 std::ifstream f("output/rank4_audit/M_integer_data.txt");int n;f>>n;if(n!=22)throw std::runtime_error("dimension");
 Mat G,W;for(auto*m:{&G,&W})for(V&r:*m)for(I&x:r)f>>x;if(!f)throw std::runtime_error("input");
 std::ifstream t("output/rank4_audit/binary_obstruction_types.txt");std::map<std::pair<I,I>,std::set<std::pair<I,I>>> wanted;I q,r,a,b;int types=0;
 while(t>>q>>r>>a>>b){wanted[{q,r}].insert({a,b});types++;}if(types==0)throw std::runtime_error("empty types");
 std::ifstream s("output/rank4_audit/M_shell12.txt");int nv;s>>nv>>n;if(n!=22||nv<=0)throw std::runtime_error("shell header");
 std::map<I,std::vector<V>> sphere;for(int j=0;j<nv;j++){V v;for(I&x:v)s>>x;I q=dot(v,mv(G,v));if(q%2||q<=0||q>24)throw std::runtime_error("shell norm");sphere[q/2].push_back(v);}if(!s)throw std::runtime_error("truncated sphere");
 std::ifstream fRep("output/rank4_audit/M_shell12_orbit_reps.txt");I q0;I comparisons=0,hits=0;int reps=0;
 while(fRep>>q0){V u;for(I&x:u)fRep>>x;if(!fRep)throw std::runtime_error("truncated representative");reps++;
  V gu=mv(G,u),gwu=mv(G,mv(W,u));if(dot(u,gu)!=2*q0)throw std::runtime_error("representative norm");
  for(const auto&kv:wanted){if(kv.first.first!=q0)continue;
   for(const V&v:sphere.at(kv.first.second)){
    I tr=dot(gu,v),trw=dot(gwu,v);if((tr-2*trw)%23)throw std::runtime_error("nonintegral Hermitian coefficient");
    I b0=(tr-2*trw)/23;if((tr-b0)%2)throw std::runtime_error("coefficient parity");I a0=(tr-b0)/2;
    for(I sign:{-1,1}){comparisons++;if(kv.second.count({sign*a0,sign*b0})){hits++;std::cerr<<"representation "<<q0<<" "<<kv.first.second<<" "<<sign*a0<<" "<<sign*b0<<"\n";}}
   }
  }
 }
 if(reps<=0)throw std::runtime_error("incomplete representative file");
 std::cout<<"{\"complete\":true,\"types\":"<<types<<",\"representatives\":"<<reps<<",\"direct_signed_comparisons\":"<<comparisons<<",\"hits\":"<<hits<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 return hits?2:0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
