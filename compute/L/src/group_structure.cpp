#define main old_minima_main
#include "../workspace/output/minima_reaudit/witnesses/minima_orbits.cpp"
#undef main
#include <set>
int main(int argc,char**argv){try{
 need(argc==5,"group_structure group.txt ambient.txt zeta23.txt output");
 std::ifstream gf(argv[1]);int n,sz;gf>>n>>sz;need(n==22&&sz==12144,"group header");std::vector<M> gs(sz);std::unordered_map<M,int,MH> ix;
 for(int i=0;i<sz;i++){for(auto&r:gs[i])for(int&v:r)gf>>v;need(ix.emplace(gs[i],i).second,"distinct group matrices");}need(bool(gf),"group read");
 std::ifstream af(argv[2]);M G{},W{};int ng;af>>n;for(auto m:{&G,&W})for(auto&r:*m)for(int&v:r)af>>v;af>>ng;std::vector<M> gen(ng),inv(ng);for(auto&m:gen)for(auto&r:m)for(int&v:r)af>>v;need(bool(af),"ambient read");
 M id{};for(int i=0;i<22;i++)id[i][i]=1;need(gs[0]==id,"identity first");
 for(int k=0;k<ng;k++){need(mul(transpose(gen[k]),mul(G,gen[k]))==G&&mul(W,gen[k])==mul(gen[k],W),"isometry generator");inv[k]=id;for(int j=0;mul(inv[k],gen[k])!=id;j++){need(j<1000,"finite generator order");inv[k]=mul(inv[k],gen[k]);}}
 std::vector<int> seen(sz);seen[0]=1;std::vector<int> todo{0};for(size_t j=0;j<todo.size();j++)for(auto&g:gen){int i=ix.at(mul(g,gs[todo[j]]));if(!seen[i]){seen[i]=1;todo.push_back(i);}}need(todo.size()==gs.size(),"full generated group");
 std::vector<int> cl(sz,-1),reps,counts;for(int i=0;i<sz;i++)if(cl[i]<0){int c=reps.size();reps.push_back(i);cl[i]=c;std::vector<int> q{i};for(size_t j=0;j<q.size();j++)for(int k=0;k<ng;k++){int a=ix.at(mul(gen[k],mul(gs[q[j]],inv[k])));if(cl[a]<0){cl[a]=c;q.push_back(a);}need(cl[a]==c,"conjugacy partition");}counts.push_back(q.size());}
 std::filesystem::create_directories(argv[4]);std::ofstream cp(std::string(argv[4])+"/conjugacy.gp");cp<<"CLASS_REPS=[";for(size_t i=0;i<reps.size();i++){if(i)cp<<',';cp<<reps[i]+1;}cp<<"];CLASS_SIZES=[";for(size_t i=0;i<counts.size();i++){if(i)cp<<',';cp<<counts[i];}cp<<"];CLASS_ID=[";for(int i=0;i<sz;i++){if(i)cp<<',';cp<<cl[i]+1;}cp<<"];\n";
 std::ifstream zf(argv[3]);zf>>n>>n;M z{};for(auto&r:z)for(int&v:r)zf>>v;need(bool(zf)&&ix.count(z),"zeta23 in group");
 auto syl=[&](M const&m){std::vector<int>s;M a=id;for(int j=1;j<=23;j++){a=mul(m,a);if(j==23)need(a==id,"order23");else{s.push_back(ix.at(a));need(a!=id,"exact order23");}}std::sort(s.begin(),s.end());return s;};
 std::vector<std::vector<int>> sy{syl(z)};std::map<std::vector<int>,int> syix;syix[sy[0]]=0;std::vector<std::vector<int>> perms(ng);
 for(size_t j=0;j<sy.size();j++)for(int k=0;k<ng;k++){auto s=syl(mul(gen[k],mul(gs[sy[j][0]],inv[k])));auto [it,b]=syix.emplace(s,sy.size());if(b)sy.push_back(s);perms[k].push_back(it->second);}
 need(sy.size()==24,"24 Sylow23 groups");std::ofstream out(std::string(argv[4])+"/sylow_action.json");out<<"{\"group_order\":12144,\"conjugacy_classes\":"<<reps.size()<<",\"sylow23_subgroups\":[";for(size_t i=0;i<sy.size();i++){if(i)out<<',';out<<'[';for(int j=0;j<22;j++){if(j)out<<',';out<<sy[i][j];}out<<']';}out<<"],\"generator_permutations\":[";for(int k=0;k<ng;k++){if(k)out<<',';out<<'[';for(int j=0;j<24;j++){if(j)out<<',';out<<perms[k][j];}out<<']';}out<<"]}\n";
 std::cout<<"GROUP_STRUCTURE_COMPLETE order="<<sz<<" classes="<<reps.size()<<" sylow23="<<sy.size()<<"\n";
 }catch(std::exception const&e){std::cerr<<e.what()<<'\n';return 1;}}
