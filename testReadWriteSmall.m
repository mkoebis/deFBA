load('modelDemo.mat')
model2 = model;

model = addMetaboliteAtPosition(model,'O2ext','O2ext',12);
model.sizeXmet = model.sizeXmet+1;
model.S(12,16) = -1;
model.genes = model.enz;
model.rxnGeneMat = model.rxnEnzRules;

% add compartments
cyt = [1:1:model.sizeXmet,model.sizeXmet+model.sizeYmet+1:1:length(model.mets)];
model.metCompartments = cell(length(model.mets),1);
for i=1:length(cyt)
    model.metCompartments{cyt(i)} = 'c';
end

ext = model.sizeXmet:1:model.sizeXmet+model.sizeYmet;
for i=1:length(ext)
    model.metCompartments{ext(i)} = 'e';
end

compL = {'c','e'};
compN = {'cytosol','extracellular'};

model.c = zeros(model.noRxn,1);
model.lb  = zeros(model.noRxn,1);
model.ub  = zeros(model.noRxn,1);

% add initial biomass
model2.ProtQuotaProdID = 'Sprod';
result = run_deFBA(model2,'cplex');
strsol = toStruct(result,model2);

model.initialBiomass = strsol.p0;

model.carbonSources = {'Carb1','Carb2','Dext','Eext','Fext','Hext'};

model.gprComp = cell(model.noRxn,1);
for i=1:model.noRxn
    model.gprComp(i) = model.enz(model.rxnEnzRules(i,:));
end


sbmlModel = writedeFBAModel(model,'modelDemo.xml',compL,compN,'e');

model = readdeFBAModel('modelDemo.xml','e');
model.ribosomeID = 'RIB';
model.ProtQuotaProdID = 'Sprod';
model.tf = 80;
model.N = 40;
r = run_deFBA(model,'cplex');
s = toStruct(r,model);

plot(s.y,'LineWidth',2)
hold on
plot(strsol.y,'--','LineWidth',2)