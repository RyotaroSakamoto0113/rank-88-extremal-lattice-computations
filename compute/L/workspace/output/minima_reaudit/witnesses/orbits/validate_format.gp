default(parisize,512000000);
read("output/minima_reaudit/witnesses/orbits/group_matrices.gp");
read("output/minima_reaudit/witnesses/orbits/orbit_representatives.gp");
print([#MINIMA_GROUP,#MINIMA_ORBIT_REPS,matsize(MINIMA_GROUP[1]),MINIMA_GROUP[1]==matid(22)]);
quit;
