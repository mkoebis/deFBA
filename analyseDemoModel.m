sol = run_deFBA(model,'cplex');
strsol = toStruct(sol,model)
strsol.y
plot(strsol.y)