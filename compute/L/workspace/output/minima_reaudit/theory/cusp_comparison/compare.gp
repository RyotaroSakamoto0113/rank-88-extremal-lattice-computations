read("/Users/ryotarosakamoto/.codex/.chatgpt-projects/g-p-6a7359f42e50819182181e61199f2684/output/minima_reaudit/runs/20260912-223806-1hi7u01i/bootstrap/cusps.gp");
BCC=CC;BGG=GG;BCI=CI;
OUT="/Users/ryotarosakamoto/.codex/.chatgpt-projects/g-p-6a7359f42e50819182181e61199f2684/output/minima_reaudit/theory/cusp_comparison";
read("/Users/ryotarosakamoto/.codex/.chatgpt-projects/g-p-6a7359f42e50819182181e61199f2684/output/minima_reaudit/runs/20260912-224127-d3-2ghz11k4/src/common.gp");
need(CC==BCC && GG==BGG && CI==BCI,"production d3/bootstrap cusp mismatch");closecert();print("D3_BOOTSTRAP_CUSPS_EXACTLY_EQUAL");quit;
