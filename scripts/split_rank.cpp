#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
using Clock=chrono::steady_clock;using Mat=vector<vector<int>>;
int prime=461,n=88;
void need(bool b,const string&s){if(!b)throw runtime_error(s);}
int mod(long long a){a%=prime;if(a<0)a+=prime;return a;}
int power(int a,int b){int r=1;for(;b;b>>=1,a=mod(a*a))if(b&1)r=mod(r*a);return r;}
Mat multiply(const Mat&a,const Mat&b){Mat c(a.size(),vector<int>(b[0].size()));for(size_t i=0;i<a.size();i++)for(size_t k=0;k<b.size();k++)if(a[i][k])for(size_t j=0;j<b[0].size();j++)c[i][j]=mod(c[i][j]+a[i][k]*b[k][j]);return c;}
int determinant(Mat a){int r=1;for(size_t k=0;k<a.size();k++){size_t i=k;while(i<a.size()&&!a[i][k])i++;if(i==a.size())return 0;if(i!=k){swap(a[i],a[k]);r=mod(-r);}r=mod(r*a[k][k]);int inv=power(a[k][k],prime-2);for(i=k+1;i<a.size();i++){int c=mod(a[i][k]*inv);for(size_t j=k;j<a.size();j++)a[i][j]=mod(a[i][j]-c*a[k][j]);}}return r;}
struct Block{
 vector<pair<int,int>> pairs;vector<vector<int>> basis;vector<int> ids;
 bool insert(vector<int>a,int id){for(int k=0;k<(int)a.size();k++)if(a[k]){if(basis[k].empty()){int inv=power(a[k],prime-2);for(int&v:a)v=mod(v*inv);basis[k]=a;ids.push_back(id);return true;}int c=a[k];for(int j=k;j<(int)a.size();j++)a[j]=mod(a[j]-c*basis[k][j]);}return false;}
};
template<class T>void jsonvec(ostream&o,const vector<T>&v){o<<"[";for(size_t i=0;i<v.size();i++){if(i)o<<",";o<<v[i];}o<<"]";}
template<class T>void jsonmat(ostream&o,const vector<vector<T>>&m){o<<"[";for(size_t i=0;i<m.size();i++){if(i)o<<",";jsonvec(o,m[i]);}o<<"]";}
int main(int argc,char**argv){try{
 need(argc>=4,"usage: split_rank integer_input candidates.txt certificate.json [prime]");if(argc>4)prime=stoi(argv[4]);need(prime>115&&(prime-1)%115==0&&prime<10000,"prime bounds");for(int d=2;d*d<=prime;d++)need(prime%d,"not prime");
 auto start=Clock::now();ifstream in(argv[1]);int nn;in>>nn;need(nn==n,"dimension");vector<vector<long long>> G(n,vector<long long>(n));Mat m5(n,vector<int>(n)),m23=m5;
 for(auto&r:G)for(auto&v:r){in>>v;need(abs(v)<1000000000,"Gram magnitude guard");}for(Mat*m:{&m5,&m23})for(auto&r:*m)for(auto&v:r){long long z;in>>z;v=mod(z);}need(bool(in),"input complete");
 Mat g=multiply(m5,m23);int zeta=2;while(!(power(zeta,115)==1&&power(zeta,5)!=1&&power(zeta,23)!=1))zeta++;need(zeta<prime,"root not found");
 vector<int>chars;for(int a=0;a<115;a++)if(gcd(a,115)==1)chars.push_back(a);need(chars.size()==88,"characters");
 Mat T;int tdet=0,tries=0;uint64_t rng=20260913;
 while(!tdet){need(++tries<100,"left cyclic covector failure");vector<int>w(n),weights(n,1),invl(n);for(int&x:w){rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;x=rng%prime;}T.assign(n,vector<int>(n));for(int a=0;a<n;a++)invl[a]=power(power(zeta,chars[a]),prime-2);
  for(int k=0;k<115;k++){for(int a=0;a<n;a++){for(int j=0;j<n;j++)T[a][j]=mod(T[a][j]+weights[a]*w[j]);weights[a]=mod(weights[a]*invl[a]);}vector<int>next(n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)next[j]=mod(next[j]+w[i]*g[i][j]);w=next;}
  tdet=determinant(T);
 }
 Mat tg=multiply(T,g);for(int a=0;a<n;a++)for(int j=0;j<n;j++)need(tg[a][j]==mod(power(zeta,chars[a])*T[a][j]),"left eigenbasis identity");
 vector<Block>blocks(115);for(int i=0;i<n;i++)for(int j=i;j<n;j++)blocks[(chars[i]+chars[j])%115].pairs.emplace_back(i,j);for(auto&b:blocks)b.basis.resize(b.pairs.size());
 auto diagend=Clock::now();ifstream seedsin(argv[2]);int amount,dimension;seedsin>>amount>>dimension;need(dimension==n&&amount>=0,"seed header");vector<vector<long long>>seeds;vector<int>sourceids;int totalrank=0,processed=0;
 for(int id=0;id<amount&&totalrank<3916;id++){
  vector<long long>x(n);for(auto&v:x){seedsin>>v;need(abs(v)<1000000000,"seed magnitude guard");}need(bool(seedsin),"seed read");__int128 norm=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++)norm+=(__int128)x[i]*G[i][j]*x[j];need(norm==8,"seed norm is not8");processed++;
  vector<int>y(n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)y[i]=mod(y[i]+T[i][j]*mod(x[j]));int proposed=seeds.size(),gain=0;
  for(auto&b:blocks){if(b.ids.size()==b.pairs.size())continue;vector<int>row;for(auto ij:b.pairs)row.push_back(mod(y[ij.first]*y[ij.second]));if(b.insert(row,proposed))gain++;}
  if(gain){seeds.push_back(x);sourceids.push_back(id);totalrank+=gain;cout<<"SEED "<<seeds.size()<<" source="<<id<<" rank="<<totalrank<<" gain="<<gain<<" invariant="<<blocks[0].ids.size()<<"\n";}
 }
 auto rankend=Clock::now();ofstream out(argv[3]);out<<"{\"status\":\""<<(totalrank==3916?"full_rank":"partial")<<"\",\"p\":"<<prime<<",\"zeta\":"<<zeta<<",\"characters\":";jsonvec(out,chars);out<<",\"T\":";jsonmat(out,T);out<<",\"seeds\":";jsonmat(out,seeds);out<<",\"source_candidate_ids\":";jsonvec(out,sourceids);out<<",\"rank\":"<<totalrank<<",\"blocks\":[";
 for(int k=0;k<115;k++){if(k)out<<",";out<<"{\"character\":"<<k<<",\"dimension\":"<<blocks[k].pairs.size()<<",\"pairs\":[";for(size_t j=0;j<blocks[k].pairs.size();j++){if(j)out<<",";out<<"["<<blocks[k].pairs[j].first<<","<<blocks[k].pairs[j].second<<"]";}out<<"],\"seed_ids\":";jsonvec(out,blocks[k].ids);out<<"}";}
 out<<"],\"timing\":{\"diagonalization_seconds\":"<<chrono::duration<double>(diagend-start).count()<<",\"rank_seconds\":"<<chrono::duration<double>(rankend-diagend).count()<<",\"processed_candidates\":"<<processed<<",\"covector_attempts\":"<<tries<<"}}\n";out.close();
 cout<<"SPLIT_RANK_RESULT "<<totalrank<<" seeds="<<seeds.size()<<" p="<<prime<<" zeta="<<zeta<<" seconds="<<chrono::duration<double>(Clock::now()-start).count()<<"\n";return totalrank==3916?0:2;
 }catch(const exception&e){cerr<<e.what()<<"\n";return 1;}}
