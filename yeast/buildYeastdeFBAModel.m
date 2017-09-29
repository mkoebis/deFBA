function deFBAmodel = buildYeastdeFBAModel
    % 	load('yeast6.06_reduced_correctNames.mat');
    load('reducedYeast6.mat');
    
    % remove steady-state total biomass and glucose, galactose, ethanol
    % exchange
    model = removeRxns(model,{'r_2111', 'r_2133','r_1710','r_1714','r_1761'});
    model = removeMetabolites(model,'s_0450[c03]');
    % remove exchanges
    model = removeRxns(model,{'r_1548','r_1654','r_1672','r_1788','r_1808',...
        'r_1832','r_1915','r_1947','r_1992','r_2005','r_2060','r_2100','r_2106','r_2134','r_2137'});
    model = removeMetabolites(model,{'s_0032[c06]','s_0420[c06]','s_0458[c06]','s_0702[c06]','s_0766[c06]',...
        's_0796[c06]','s_1061[c06]','s_1154[c06]','s_1277[c06]','s_1324[c06]',...
        's_1468[c06]','s_0805[c06]','s_1571[c06]','s_2766[c01]','s_2768[c01]'});
    
%     % remove dead ends
%     model = removeMetabolites(model,{'s_0164[c11]','s_0383[c03]','s_0685[c03]','s_0868[c07]','s_0889[c07]',...
%         's_1096[c03]','s_1117[c07]','s_1126[c07]','s_1129[c07]','s_1132[c07]','s_1135[c07]','s_1138[c07]','s_1141[c07]','s_1144[c07]','s_1525[c09]'});
%     % remove corresp. reactions
%     model = removeRxns(model,{'r_1574','r_0139','r_0369','r_0580','r_0585','r_1509','r_1510','r_1511','r_1512','r_1513','r_1516','r_1517','r_2081','r_2108'});

    % remove putrescine cycling
    model = removeRxns(model,{'r_1250', 'r_1251'});
    % dUTP diphosphatase should only go forward
    model.rev(findRxnIDs(model,'r_0364')) = 0;
    % http://www.yeastgenome.org/reference/S000086806/overview
    model.rev(findRxnIDs(model,'r_1130')) = 0;
    model.rev(findRxnIDs(model,'r_1131')) = 0;
    % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3026612/#R53
    model.rev(findRxnIDs(model,'r_0735')) = 0;
    model.rev(findRxnIDs(model,'r_0736')) = 0;
    % https://www.ncbi.nlm.nih.gov/pubmed/10979554
    model.rev(findRxnIDs(model,'r_0143')) = 0;
    % http://onlinelibrary.wiley.com/doi/10.1046/j.1365-2958.2003.03810.x/full
    model.rev(findRxnIDs(model,'r_1245')) = 0;
    % pyruvate kinase is irreversible
    model.rev(findRxnIDs(model,'r_0962')) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% yeast 7 changes:
    %     % remove H+ in glycine-cleavage complex rxn [r_0507] for charge balancing
    model.S(findMetIDs(model,'s_0799[c11]'),findRxnIDs(model,'r_0507')) = 0;

    % GDH2 is irreversible in the backwards direction
    model.S(:,findRxnIDs(model,'r_0470')) = -model.S(:,findRxnIDs(model,'r_0470'));
    model.rev(findRxnIDs(model,'r_0470')) = 0;
    % r_1110 does not cotransport protons
    model.S(findMetIDs(model,{'s_0794[c03]','s_0799[c11]'}),findRxnIDs(model,'r_1110')) = 0;
    % ATP synthase r_0226 moves 4 cytoplasmic protons
    model.S(findMetIDs(model,{'s_0794[c03]','s_0799[c11]'}),findRxnIDs(model,'r_0226')) = [-4,3];
    % NADH oxidoreductase r_0773 is not proton translocating
    model.S(findMetIDs(model,{'s_0794[c03]','s_0799[c11]'}),findRxnIDs(model,'r_0773')) = [0,1];
    % r_0718 is NADP-dependent, not NAD
    model.S(findMetIDs(model,{'s_1200[c11]','s_1205[c11]'}),findRxnIDs(model,'r_0718')) = 0;
    model.S(findMetIDs(model,{'s_1210[c11]','s_1214[c11]'}),findRxnIDs(model,'r_0718')) = [-1,1];
    % proton stoichiometries for r_0438 and r_0439
    model.S(findMetIDs(model,{'s_0794[c03]','s_0799[c11]'}),findRxnIDs(model,'r_0438')) = [4,-8];
    model.S(findMetIDs(model,{'s_0794[c03]','s_0799[c11]'}),findRxnIDs(model,'r_0439')) = [4,-2];
    % r_0678 requires ATP
    model.S(findMetIDs(model,{'s_0434[c03]','s_0803[c03]','s_0423[c03]','s_0633[c03]'}),findRxnIDs(model,'r_0678'))= [-1,0,1,1];
    % r_1099 change metabolite from malate to 2-oxoglutarate
    model.S(findMetIDs(model,{'s_0066[c03]','s_0180[c03]','s_0068[c11]','s_0182[c11]'}),findRxnIDs(model,'r_1099'))= [0,-1,0,1];

    % problem: ADH2 is not retained in the reduced model, which makes it
    % that ethanol can't be reutilized
    % => we have to add it
    model = addReaction(model,{'r_0163','alcohol dehydrogenase (ethanol to acetaldehyde)'},...
        {'s_0680[c03]','s_1198[c03]','s_0359[c03]','s_0794[c03]','s_1203[c03]'},...
        [-1,-1,1,1,1],false,'YMR303C');

    % fructose bisphosphatase is also not retained, but is a key enzyme of
    % gluconeogenesis => add it back
    model = addReaction(model,{'r_0449','fructose-bisphosphatase'},...
        {'s_0434[c03]','s_1271[c03]','s_0394[c03]','s_0456[c03]','s_1360[c03]'},...
        [-1,-1,1,1,1],false,'YLR377C');

    % pep carboxykinase is also not retained, but is a key enzyme of
    % gluconeogenesis => add it back
    model = addReaction(model,{'r_0884','phosphoenolpyruvate carboxykinase'},...
        {'s_0555[c03]','s_0803[c03]','s_0557[c03]','s_1322[c03]'},...
        [-1,-1,1,1],false,'YKR097W');
    % allow growth on ethanol
    model = addMetabolite(model,'s_0361[c11]','acetaldehyde [mitochondrion]');
    model = addReaction(model,{'r_1632','acetaldehyde transport'},...
        {'s_0361[c11]','s_0359[c03]'},...
        [-1,1],true,'');

    model = addMetabolite(model,'s_0365[c11]','acetate [mitochondrion]');
    model = addReaction(model,{'r_0174','aldehyde dehydrogenase (acetylaldehyde, NAD)'},...
        {'s_0361[c11]','s_0807[c11]','s_1200[c11]','s_0365[c11]','s_0799[c11]','s_1205[c11]'},...
        [-1,-1,-1,1,2,1],false,'YOR374W');

    model = addReaction(model,{'r_0175','aldehyde dehydrogenase (acetylaldehyde, NADP)'},...
        {'s_0361[c11]','s_0807[c11]','s_1210[c11]','s_0365[c11]','s_0799[c11]','s_1214[c11]'},...
        [-1,-1,-1,1,2,1],false,'YER073W OR YOR374W');

    model = addMetabolite(model,'s_0424[c11]','AMP [mitochondrion]');
    model = addReaction(model,{'r_0113','acetyl-CoA synthetase'},...
        {'s_0365[c11]','s_0437[c11]','s_0532[c11]','s_0376[c11]','s_0424[c11]','s_0636[c11]'},...
        [-1,-1,-1,1,1,1],false,'YAL054C');

    model = addReaction(model,{'r_0149','adenylate kinase'},...
        {'s_0424[c11]','s_0437[c11]','s_0397[c11]'},...
        [-1,-1,2],false,'YER170W');

    % add back reactions that allow storage usage
    model = addReaction(model,{'r_0463','glucan 1,4-alpha-glucosidase'},...
        {'s_0773[c03]','s_0803[c03]','s_0563[c03]'},...
        [-1,-1,1], false,'YPR184W');
    model = addReaction(model,{'r_0511','glycogen phosphorylase'},...
        {'s_0773[c03]','s_1322[c03]','s_0567[c03]'},...
        [-1,-1,1],false,'YPR160W');
    model = addReaction(model,{'r_0194','alpha,alpha-trehalase'},...
        {'s_0803[c03]','s_1520[c03]','s_0563[c03]'},...
        [-1,-1,2], false,'YDR001C');
    
    % add back cytosolic glycerol 3-phosphate dehydrogenase
    model = addReaction(model,{'r_0491','glycerol-3-phosphate dehydrogenase (NAD)'},...
        {'s_0629[c03]','s_0794[c03]','s_1203[c03]','s_0767[c03]','s_1198[c03]'},...
        [-1,-1,-1,1,1],false,'YOL059W OR YDL022W');


    % remove gene association for water uptake: uniprot says it's
    % unfunctional in most strains, so reaction should be spontaneous
    model = changeGeneAssociation(model,'r_1277','');

    % get reversibility from lb
    model.rev(model.lb==0)=0;

    % reversibility settings from yeast7:
    model.rev(findRxnIDs(model,{'r_1811','r_2045'})) = 1;
    model.S(:,findRxnIDs(model,'r_0217')) = -model.S(:,findRxnIDs(model,'r_0217'));
    model.rev(findRxnIDs(model,'r_0217')) = 0;
    
    % add maintenance reaction
    model = addMaintenanceYeast(model);

    % remove ors and duplicate reactions instead
    model = duplicateOrRxns(model);

    % exchange should look in the same direction
    model.S(:,findRxnIDs(model,{'r_1697'})) = ...
        - model.S(:,findRxnIDs(model,{'r_1697'}));

    % add quota reactions; growth-associated maintenance included
    % proteins
    % energy required for polymerization: 16.965 mmol ATP/g cell (Verduyn,
    % 1991, Jochen Foerster et al 2003 supplement)
    deFBAmodel = addMetabolite(model,'B1[c03]','Protein quota');
    deFBAmodel = addReaction(deFBAmodel,{'synth_BM0001','Production of protein quota'},...
        {'s_0404[c03]','s_0428[c03]','s_0430[c03]','s_0432[c03]',...
        's_0542[c03]','s_0747[c03]','s_0748[c03]','s_0757[c03]',...
        's_0832[c03]','s_0847[c03]','s_1077[c03]','s_1099[c03]',...
        's_1148[c03]','s_1314[c03]','s_1379[c03]','s_1428[c03]',...
        's_1491[c03]','s_1527[c03]','s_1533[c03]','s_1561[c03]',...
        'B1[c03]','s_0434[c03]','s_0803[c03]',...
        's_1582[c03]','s_1583[c03]','s_1585[c03]','s_1587[c03]',...
        's_1589[c03]','s_1590[c03]','s_1591[c03]','s_1593[c03]',...
        's_1594[c03]','s_1596[c03]','s_1598[c03]','s_1600[c03]',...
        's_1602[c03]','s_1604[c03]','s_1606[c03]','s_1607[c03]',...
        's_1608[c03]','s_1610[c03]','s_1612[c03]','s_1614[c03]',...
        's_0394[c03]','s_1322[c03]','s_0794[c03]'},...
        [-0.9839201541, -0.3446294001, -0.2181008711, -0.6380040232,...
        -0.0141540388, -0.2260357111, -0.6472255939, -0.6227777087,...
        -0.1421837537, -0.4132550429, -0.6356450167, -0.6137705931,...
        -0.1087287529, -0.2871554242, -0.3532076054, -0.3975998181,...
        -0.4104671262, -0.060905258, -0.2187442365, -0.5674482841,...
        1, -36.3823134562, -36.3823134562,...
        0.9839201541, 0.3446294001, 0.2181008711, 0.6380040232,...
        0.0141540388, 0.2260357111, 0.6472255939, 0.6227777087,...
        0.1421837537, 0.4132550429, 0.6356450167, 0.6137705931,...
        0.1087287529, 0.2871554242, 0.3532076054, 0.3975998181,...
        0.4104671262, 0.060905258, 0.2187442365, 0.5674482841,...
        36.3823134562,36.3823134562,36.3823134562], false);

    %RNA
    % energy required for polymerization: 1.638 mmol ATP/g cell (Verduyn,
    % 1991, Jochen Foerster et al 2003 supplement)
    deFBAmodel = addMetabolite(deFBAmodel,'B2[c03]','RNA quota');
    deFBAmodel = addReaction(deFBAmodel,{'synth_BM0002','Production of RNA quota'},...
        {'s_0423[c03]','s_0526[c03]','s_0782[c03]','s_1545[c03]',...
        'B2[c03]','s_0434[c03]','s_0803[c03]',...
        's_0394[c03]','s_1322[c03]','s_0794[c03]'},...
        [-0.6912546756, -0.6717192174, -0.6912546756, -0.9001338059,...
        1, -24.6146773624, -24.6146773624,...
        24.6146773624,24.6146773624,24.6146773624], false);

    %DNA
    % energy required for polymerization: 0.104 mmol ATP/g cell (Verduyn,
    % 1991, Jochen Foerster et al 2003 supplement)
    deFBAmodel = addMetabolite(deFBAmodel,'B3[c03]','DNA quota');
    deFBAmodel = addReaction(deFBAmodel,{'synth_BM0003','Production of DNA quota'},...
        {'s_0584[c03]','s_0589[c03]','s_0615[c03]','s_0649[c03]',...
        'B3[c03]','s_0434[c03]','s_0803[c03]',...
        's_0394[c03]','s_1322[c03]','s_0794[c03]'},...
        [-0.9193796026, -0.6129197351, -0.6129197351, -0.9193796026,...
        1, -26.5598551875,-26.5598551875,...
        26.5598551875,26.5598551875,26.5598551875], false);

    % Cell wall
    % energy required for polymerization: 5.2096 mmol ATP/g cell (Verduyn,
    % 1991, Jochen Foerster et al 2003 supplement)
    deFBAmodel = addMetabolite(deFBAmodel,'B4[c03]','Cell wall quota');
    deFBAmodel = addReaction(deFBAmodel,{'synth_BM0004','Production of cell wall quota'},...
        {'s_1107[c03]','s_0002[c03]',...
        'B4[c03]','s_0434[c03]','s_0803[c03]',...
        's_0394[c03]','s_1322[c03]','s_0794[c03]'},...
        [-2.5670649287, -3.6057745774,...
        1, -16.5532633404,-16.5532633404,...
        16.5532633404, 16.5532633404, 16.5532633404], false);

    % Lipids
    deFBAmodel = addMetabolite(deFBAmodel,'B5[c03]','Lipids quota');
    deFBAmodel = addReaction(deFBAmodel,{'synth_BM0005','Production of lipids quota'},...
        {'s_0089[c03]','s_0122[c03]','s_0535[c03]','s_0657[c03]',...
        's_0662[c03]','s_0666[c03]','s_0672[c03]','s_0694[c03]',...
        's_0700[c03]','s_1059[c03]','s_1337[c03]','s_1346[c03]',...
        's_1351[c03]','s_1524[c03]','s_1569[c03]',...
        'B5[c03]'},...
        [-0.1852844279, -0.0067816523, -0.0504990892, -0.0116256896,...
        -0.0151376167, -0.6781652264, -0.0983339578, -0.0249467923,...
        -0.0138055064, -0.0038752299, -0.0451706481, -0.3487706879,...
        -0.0844073505, -0.0945798289, -0.001816514,...
        1], false);

    % Small molecules
    deFBAmodel = addMetabolite(deFBAmodel,'B6[c03]','Small molecules quota');
    deFBAmodel = addReaction(deFBAmodel,{'synth_BM0006','Production of small molecules quota'},...
        {'s_1405[c03]','s_1467[c03]',...
        'B6[c03]'},...
        [-0.4315988987, -8.71916967, 1], false);

       
    deFBAmodel = orderMetsRxnsYeast(deFBAmodel);

    % translation rate: http://www.ncbi.nlm.nih.gov/pubmed/1097403%2C1089627?dopt=Abstract
    deFBAmodel.ribosomeRate = 10;

    deFBAmodel = formulateBiomassYeast(deFBAmodel,'AminoAcidCountsAll.csv');
    
    deFBAmodel.quotaInitial = [0.466298; 0.066545662; 0.003915684;...
        0.3147174; 0.0082575747; 0.0022937964; zeros(length(deFBAmodel.enz),1)];
    % correct steady-state protein biomass to only 0.5947 (only non-metabolic proteins, since the rest are modeled explicitly)
    % number from data in http://www.ncbi.nlm.nih.gov/pubmed/18820680
    deFBAmodel.quotaInitial(1) = deFBAmodel.quotaInitial(1)*0.5947;

    % remove unnecessary fields
    % 	deFBAmodel = rmfield(deFBAmodel,'c');
    deFBAmodel = rmfield(deFBAmodel,'b');
    deFBAmodel = rmfield(deFBAmodel,'description');


    deFBAmodel.noRxn = length(deFBAmodel.rxns);

    for i=1:length(deFBAmodel.rxns)
        assert(length(find(ismember(deFBAmodel.rxns,deFBAmodel.rxns(i))))==1)
    end
    % set Kcat_f and Kcat_r using rxnEnzRules and reversibility for the moment
    deFBAmodel.Kcat_f = deFBAmodel.rxnEnzRules*700000; % ribosome Kcat set separately
    % 	deFBAmodel.Kcat_r = deFBAmodel.rxnEnzRules*700000;

    % load EC numbers
    T = readtable('ECnumbers_yeast_reduced.csv');
    deFBAmodel.rxnECNumbers = cell(deFBAmodel.noRxn,1);
    for i=1:length(T.rxn)
        idx = findRxnIDs(deFBAmodel,T.rxn(i));
        if idx~=0
            if ~isempty(T.EC{i}) && isempty(strfind(T.EC{i},'.-'))
                deFBAmodel.rxnECNumbers{idx} = T.EC{i};
            else
                deFBAmodel.rxnECNumbers{idx} = '';
            end
        end
    end
    for i=1:deFBAmodel.noRxn
        if ~ischar(deFBAmodel.rxnECNumbers{i})
            deFBAmodel.rxnECNumbers{i} = '';
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5got up to here
    % assign Kcat from database
    T = readtable('data_kcat_yeast_phylogeny_YeastPreferred.csv');
    deFBAmodel = assignDataKcatYeast(deFBAmodel,T.EC,T.kbest);

    % set Kcat of ribosome
    en = find(deFBAmodel.rxnEnzRules(:,end));
    deFBAmodel.Kcat_f(en,end) = deFBAmodel.ribosomeRate*3600./(deFBAmodel.rxnEnzRules(en,end).*deFBAmodel.enzLength);

    % set kcat of glucose and galactose transporters
    transp = {'r_1166','r_1135'};
    for j=1:length(transp)
        idx = find(strncmp(deFBAmodel.rxns,transp{j},6));
        for i=1:length(idx)
            deFBAmodel.Kcat_f(findRxnIDs(deFBAmodel,deFBAmodel.rxns(idx(i))),deFBAmodel.Kcat_f(findRxnIDs(deFBAmodel,deFBAmodel.rxns(idx(i))),:)~=0) = 200*3600;
        end
    end
    % correct for nonmetabolic protein quota reaction
    AAnames = {'s_0404[c03]','s_0428[c03]','s_0430[c03]','s_0432[c03]',...
        's_0542[c03]','s_0747[c03]','s_0748[c03]','s_0757[c03]',...
        's_0832[c03]','s_0847[c03]','s_1077[c03]','s_1099[c03]',...
        's_1148[c03]','s_1314[c03]','s_1379[c03]','s_1428[c03]',...
        's_1491[c03]','s_1527[c03]','s_1533[c03]','s_1561[c03]'};
    AAids = findMetIDs(deFBAmodel,AAnames);
    idxBM1 = findRxnIDs(deFBAmodel,'synth_BM0001');
    deFBAmodel.spontaneousRxn(idxBM1) = 0;
    deFBAmodel.rxnEnzRules(idxBM1,end) = 1;
    deFBAmodel.rxnGeneMat(idxBM1,:) = deFBAmodel.rxnGeneMat(end,:);% gene association for protein quota is the same as for ribosome
    deFBAmodel.Kcat_f(idxBM1,end) = deFBAmodel.ribosomeRate*3600/sum(-deFBAmodel.S(AAids,idxBM1));

    % set N,K,tf,epsilon,beta
    % 	deFBAmodel.tf = 3;
    deFBAmodel.storageWeight = [180.1559,342.296]*1e-3;
    deFBAmodel.N = 10;
    deFBAmodel.K = 1;
    deFBAmodel.epsilon = 1e-4;
    deFBAmodel.beta = 0;
    deFBAmodel.units.quotaInitial = 'grams/gDW_cell';
    deFBAmodel.units.S = 'for usual metabolites: mmol/gDW; for proteins: mmol/gDW_cell; for quota: g/gDW';
    deFBAmodel.units.lbub = 'for usual rxns: mmol/(gDW_cell*h); for BM0001-BM0007: 1/(gDW_cell*h)';
    deFBAmodel.units.Y0 = 'mmol/gDW_cell';
    deFBAmodel.units.Kcat_fKcat_r = '1/h';
    deFBAmodel.units.tf = 'h';
    deFBAmodel.units.N = 'number of discretization points';
    deFBAmodel.units.beta = '1/h';
    deFBAmodel.units.dryWeight = 'g';
    deFBAmodel.rxnExtraBounds = {};
    deFBAmodel.maintenanceID = 'ATP_maintenance';
    deFBAmodel.maintenanceValue = 0;
    deFBAmodel.name = 'Yeast_deFBA_model';
    deFBAmodel.note = 'This model uses the ram standard 1.0';
    assert(isempty(find(deFBAmodel.Kcat_f<0)))
    deFBAmodel.tf = 10;
    deFBAmodel.Y0 = [0 0 1000 0 0];
    deFBAmodel.quotaWeights = ones(deFBAmodel.sizeQuotaMet,1);
    deFBAmodel.objectiveWeights = [zeros(deFBAmodel.sizeQuotaMet,1); deFBAmodel.proteinWeights];
    
    % clean up names
    for i=1:length(deFBAmodel.mets)
        deFBAmodel.mets{i} = strrep(deFBAmodel.mets{i},'_[','_');
        deFBAmodel.mets{i} = strrep(deFBAmodel.mets{i},']','');
    end
    
    for i=1:length(deFBAmodel.mets)
        deFBAmodel.mets{i} = strrep(deFBAmodel.mets{i},'[','_');
    end
    
    for i=1:length(deFBAmodel.rxns)
        deFBAmodel.rxns{i} = strrep(deFBAmodel.rxns{i},'_[','_');
        deFBAmodel.rxns{i} = strrep(deFBAmodel.rxns{i},']','');
    end
    
    for i=1:length(deFBAmodel.rxns)
        deFBAmodel.rxns{i} = strrep(deFBAmodel.rxns{i},'[','_');
    end
    
    for i=1:length(deFBAmodel.enz)
        deFBAmodel.enz{i} = strrep(deFBAmodel.enz{i},'_[','_');
        deFBAmodel.enz{i} = strrep(deFBAmodel.enz{i},']','');
    end
    
    for i=1:length(deFBAmodel.enz)
        deFBAmodel.enz{i} = strrep(deFBAmodel.enz{i},'[','_');
    end
    
    deFBAmodel.ribosomeID = deFBAmodel.mets{end};
    deFBAmodel.ProtQuotaProdID = 'synth_BM0001';
    deFBAmodel = rmfield(deFBAmodel,'lb');
    deFBAmodel = rmfield(deFBAmodel,'ub');
    deFBAmodel = rmfield(deFBAmodel,'c');
end

