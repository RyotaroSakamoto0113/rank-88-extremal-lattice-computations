// Independent exhaustive Gram enumeration over Z[(1+sqrt(-23))/2].
// Every arithmetic operation affecting inclusion/exclusion is integral.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
using I=long long;
struct F{I a=0,b=0;};
F operator+(F x,F y){return{x.a+y.a,x.b+y.b};}
F operator-(F x,F y){return{x.a-y.a,x.b-y.b};}
F operator-(F x){return{-x.a,-x.b};}
F operator*(F x,F y){return{x.a*y.a-6*x.b*y.b,x.a*y.b+x.b*y.a+x.b*y.b};}
F operator*(I n,F x){return{n*x.a,n*x.b};}
bool operator==(F x,F y){return x.a==y.a&&x.b==y.b;}
F bar(F x){return{x.a+x.b,-x.b};}
I norm(F x){return x.a*x.a+x.a*x.b+6*x.b*x.b;}
I trace(F x){return 2*x.a+x.b;}
I mod(I a,I m){I r=a%m;return r<0?r+m:r;}
I residue(F x){return mod(x.a+12*x.b,23);}
using Row=std::array<F,4>;
using Mat=std::array<Row,4>;
using R2=std::array<F,2>;
using M2=std::array<R2,2>;
F dot(R2 x,R2 y){return x[0]*y[0]+x[1]*y[1];}
R2 rowmul(R2 x,M2 H){return {x[0]*H[0][0]+x[1]*H[1][0],x[0]*H[0][1]+x[1]*H[1][1]};}
F hgram(R2 x,M2 H,R2 y){return dot(rowmul(x,H),{bar(y[0]),bar(y[1])});}
F hgram(Row x,Mat H,Row y){
 F out;
 for(int i=0;i<4;i++)for(int j=0;j<4;j++)out=out+x[i]*H[i][j]*bar(y[j]);
 return out;
}
Mat transform(Mat R,Mat H){
 Mat A{};
 for(int i=0;i<4;i++)for(int j=i;j<4;j++){A[i][j]=hgram(R[i],H,R[j]);A[j][i]=bar(A[i][j]);}
 return A;
}
I det3(Mat H,int a,int b,int c){
 return H[a][a].a*H[b][b].a*H[c][c].a-H[a][a].a*norm(H[b][c])-H[b][b].a*norm(H[a][c])-H[c][c].a*norm(H[a][b])+trace(H[a][b]*H[b][c]*H[c][a]);
}
std::map<I,std::vector<F>> normcache;
const std::vector<F>& normvals(I B){
 auto it=normcache.find(B);if(it!=normcache.end())return it->second;
 std::vector<F> v;
 // N(a+bw)=(a+b/2)^2+23b^2/4. These deliberately loose
 // integral boxes contain the full ellipse; final test is exact.
 I lim=0;while(23*lim*lim<=4*B)lim++;
 I al=0;while(al*al<=B)al++;
 for(I b=-lim;b<=lim;b++)for(I a=-al-std::abs(b);a<=al+std::abs(b);a++)if(norm({a,b})<=B)v.push_back({a,b});
 return normcache.emplace(B,std::move(v)).first->second;
}
std::map<std::pair<I,I>,std::vector<F>> rescache;
const std::vector<F>& resvals(I B,I res){
 auto key=std::make_pair(B,res);auto it=rescache.find(key);if(it!=rescache.end())return it->second;
 std::vector<F> v;for(F x:normvals(B))if(residue(x)==res)v.push_back(x);
 return rescache.emplace(key,std::move(v)).first->second;
}
std::map<std::array<I,3>,std::vector<F>> kcache;
const std::vector<F>& kvals(I B,F r){
 std::array<I,3>key{B,mod(r.a,23),mod(r.b,23)};
 auto it=kcache.find(key);if(it!=kcache.end())return it->second;
 std::vector<F>v;for(F x:normvals(B))if(mod(x.a-r.a,23)==0&&mod(x.b-r.b,23)==0)v.push_back(x);
 return kcache.emplace(key,std::move(v)).first->second;
}
struct Ext{R2 r;I ell;};
std::vector<Ext> extensions(M2 K,R2 C,I n,I delta){
 R2 r0=rowmul(C,K);I b=K[0][0].a,a=K[1][1].a;
 std::vector<Ext> out;
 for(F u:resvals(n*b-16,residue(r0[0])))for(F v:resvals(n*a-16,residue(r0[1]))){
  I value=a*norm(u)+b*norm(v)-trace(u*K[0][1]*bar(v)),ell=n*delta-value;
  if(ell>=552)out.push_back({{u,v},ell});
 }
 return out;
}
using Lin=std::array<I,16>;
std::vector<Mat> charts;
std::vector<Lin> direction;
Lin normform(Row v){
 Lin L{};for(int i=0;i<4;i++)L[i]=norm(v[i]);int p=4;
 for(int i=0;i<4;i++)for(int j=i+1;j<4;j++){F r=v[i]*bar(v[j]);L[p++]=trace(r);L[p++]=trace(r*F{0,1});}
 return L;
}
std::array<I,16> coords(Mat H){
 std::array<I,16>v{};for(int i=0;i<4;i++)v[i]=H[i][i].a;int p=4;
 for(int i=0;i<4;i++)for(int j=i+1;j<4;j++){v[p++]=H[i][j].a;v[p++]=H[i][j].b;}return v;
}
bool chartnorms(Mat H,I s,I T){
 auto v=coords(H);std::array<I,120>ns{};
 for(int i=0;i<120;i++){
  I n=std::inner_product(direction[i].begin(),direction[i].end(),v.begin(),I(0));ns[i]=n;
  if(n%529||n<6*529||n>(T-18)*529||n==7*529)return false;
 }
 for(int i=0;i<30;i++)if(std::min(ns[4*i]+ns[4*i+1],ns[4*i+2]+ns[4*i+3])<s*529)return false;
 return true;
}
bool chartminors(Mat H){
 for(Mat R:charts){Mat A=transform(R,H);
  for(int i=0;i<4;i++)for(int j=i+1;j<4;j++)if(A[i][i].a*A[j][j].a-norm(A[i][j])<16*529LL*529)return false;
  for(int i=0;i<4;i++)for(int j=i+1;j<4;j++)for(int k=j+1;k<4;k++)if(det3(A,i,j,k)<552*529LL*529*529)return false;
 }
 return true;
}
Mat inverse_scaled(){
 F d{-1,2};Mat R{};
 R[0]={-d,F{},-4*d,(-d)*F{1,-1}};
 R[1]={F{},-d,(-d)*F{0,-1},-4*d};
 R[2][2]={23,0};R[3][3]={23,0};return R;
}
void printF(std::ostream&o,F x){o<<"["<<x.a<<","<<x.b<<"]";}
void printMat(std::ostream&o,Mat H){
 o<<"[";for(int i=0;i<4;i++){if(i)o<<",";o<<"[";for(int j=0;j<4;j++){if(j)o<<",";printF(o,H[i][j]);}o<<"]";}o<<"]";
}
int main(int argc,char**argv){
 I h=argc>1?std::stoll(argv[1]):2,T=23*h;
 assert(h==2||h==3);
 std::string stem=argc>2?argv[2]:"output/rank4_audit/independent_h"+std::to_string(h);
 std::ifstream ci("output/rank4_audit/charts_scaled23.json");assert(ci);
 std::string raw((std::istreambuf_iterator<char>(ci)),{});for(char&c:raw)if(!(c=='-'||(c>='0'&&c<='9')))c=' ';
 std::istringstream is(raw);std::vector<I>nums;I x;while(is>>x)nums.push_back(x);assert(nums.size()==960);
 for(int k=0;k<30;k++){Mat R{};for(int i=0;i<4;i++)for(int j=0;j<4;j++)R[i][j]={nums[32*k+i+4*j],nums[32*k+i+4*j+16]};charts.push_back(R);for(Row r:R)direction.push_back(normform(r));}
 R2 C1{{{-4,0},{-1,1}}},C2{{{0,1},{-4,0}}};
 std::map<std::string,I>cnt{{"cores",0},{"branches",0},{"extension_branches",0},{"extension_pairs",0},{"k_candidates",0},{"positive",0},{"initial_minors",0},{"directions_blocks",0},{"charts",0}};
 auto start=std::chrono::steady_clock::now();
 std::ofstream targets(stem+"_targets.jsonl"),ys(stem+"_HY.jsonl");
 Mat E=inverse_scaled();
 for(I a=6;a<=T/4;a++){if(a==7)continue;
  for(I b=a;b<=T/2-a;b++){if(b==7)continue;I s=a+b;
   for(F c:normvals(a*b-16)){
    cnt["cores"]++;I delta=a*b-norm(c);if(delta*(T-s)<1104)continue;
    M2 K{{{{{b,0},c}},{{bar(c),{a,0}}}}};
    F q1=hgram(C1,K,C1),q2=hgram(C2,K,C2);assert(q1.b==0&&q2.b==0);
    for(I n1=6;n1<=T-s-6;n1++){I n2=T-s-n1;if(n1==7||n2==7||mod(n1-q1.a,23)||mod(n2-q2.a,23))continue;
     cnt["branches"]++;auto e1=extensions(K,C1,n1,delta),e2=extensions(K,C2,n2,delta);if(e1.empty()||e2.empty())continue;
     cnt["extension_branches"]++;
     for(Ext u:e1)for(Ext v:e2){
      cnt["extension_pairs"]++;
      F k0=dot(u.r,{bar(C2[0]),bar(C2[1])})+dot(C1,{bar(v.r[0]),bar(v.r[1])})-hgram(C1,K,C2);
      M2 adj{{{{{a,0},-c}},{{-bar(c),{b,0}}}}};
      F p12=hgram(u.r,adj,v.r);
      for(F k:kvals(n1*n2-16,k0)){
       cnt["k_candidates"]++;I num=u.ell*v.ell-norm(delta*k-p12);if(num<=0)continue;assert(num%delta==0);I detY=num/delta;assert(detY%529==0);I detX=detY/529;assert(detX<=529*h*h*h*h/256);
       cnt["positive"]++;
       Mat H{{{{{n1,0},k,u.r[0],u.r[1]}},{{bar(k),{n2,0},v.r[0],v.r[1]}},{{bar(u.r[0]),bar(v.r[0]),{b,0},c}},{{bar(u.r[1]),bar(v.r[1]),bar(c),{a,0}}}}};
       bool pass=true;
       for(int i=0;i<4;i++)for(int j=i+1;j<4;j++)if(H[i][i].a*H[j][j].a-norm(H[i][j])<16)pass=false;
       for(int i=0;i<4;i++)for(int j=i+1;j<4;j++)for(int l=j+1;l<4;l++)if(det3(H,i,j,l)<552)pass=false;
       if(!pass)continue;cnt["initial_minors"]++;
       if(!chartnorms(H,s,T))continue;cnt["directions_blocks"]++;
       if(!chartminors(H))continue;
       I id=cnt["charts"]++;Mat HX=transform(E,H);
       for(Row&r:HX)for(F&v:r){assert(v.a%529==0&&v.b%529==0);v.a/=529;v.b/=529;}
       targets<<"{\"id\":"<<id<<",\"det\":"<<detX<<",\"H\":";printMat(targets,HX);targets<<"}\n";
       ys<<"{\"id\":"<<id<<",\"det\":"<<detX<<",\"H\":";printMat(ys,H);ys<<"}\n";
      }
     }
    }
   }
  }
  std::cerr<<"a="<<a<<" pairs="<<cnt["extension_pairs"]<<" candidates="<<cnt["charts"]<<"\n";
 }
 double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 std::ofstream summary(stem+"_abstract_cpp.json");
 summary<<"{\"h\":"<<h<<",\"seconds\":"<<sec<<",\"counts\":{";bool first=true;for(auto kv:cnt){if(!first)summary<<",";first=false;summary<<"\""<<kv.first<<"\":"<<kv.second;}summary<<"}}\n";
 std::cout<<"COMPLETE h="<<h<<" candidates="<<cnt["charts"]<<" seconds="<<sec<<"\n";
}
