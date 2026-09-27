default(nbthreads,1);
default(parisizemax,1000000000);
INPUTS=Str(ROOT,"/workspace/output/d4_projective_J/inputs");
OUT_JSON=Str(OUT,"/verification.json");
read(Str(ROOT,"/workspace/output/d4_projective_J/src/check_initial.gp"));
