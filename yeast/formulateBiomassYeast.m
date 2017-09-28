function newModel = formulateBiomassYeast(model,aaFilename)
    %% get amino acid counts for each gene
    M = importdata(aaFilename,',',1);
    geneNames = M.textdata(2:end,1);
    AAcounts = M.data;
    % of these consider only metabolic genes
    idx = ismember(geneNames,model.genes);
    %geneNames = geneNames(idx);
    AAcounts = AAcounts(idx,:);
    geneNames = geneNames(idx);
    
    %% additionally get stoichiometry of gene products in enzymes
    T = readtable('enzymeStoichiometry.csv','ReadVariableNames',false);
    T = T(:,1:4);
    T.Properties.VariableNames = {'rxns','rxnNames','geneAssoc','stoichiometry'};
    
    %% get ribosome composition
    Tgenes = readtable('yeastRibosomalGenes.csv','ReadVariableNames',false);
    Tgenes = table2cell(Tgenes);
%     model = addGenes(model,Tgenes);
%     geneAssocRib = strjoin(Tgenes,' AND ');
    geneAssocRib = 'YGR214W AND YLR048W AND YGL123W AND YNL178W AND YLR441C AND YML063W AND YJR145C AND YHR203C AND YJR123W AND YPL090C AND YBR181C AND YOR096W AND YNL096C AND YBL072C AND YER102W AND YPL081W AND YBR189W AND YOR293W AND YMR230W AND YDR025W AND YBR048W AND YOR369C AND YDR064W AND YCR031C AND YJL191W AND YOL040C AND YJL190C AND YLR367W AND YMR143W AND YDL083C AND YML024W AND YDR447C AND YDR450W AND YML026C AND YOL121C AND YNL302C AND YHL015W AND YKR057W AND YJL136C AND YGR118W AND YPR132W AND YER074W AND YIL069C AND YGR027C AND YLR333C AND YGL189C AND YER131W AND YKL156W AND YHR021C AND YLR167W AND YOR167C AND YLR264W AND YLR388W AND YDL061C AND YLR287C-A AND YOR182C AND YOR063W AND YBR031W AND YDR012W AND YPL131W AND YML073C AND YLR448W AND YGL076C AND YPL198W AND YHL033C AND YLL045C AND YFR031C-A AND YIL018W AND YGL147C AND YNL067W AND YLR075W AND YPL220W AND YGL135W AND YPR102C AND YGR085C AND YEL054C AND YDR418W AND YDL082W AND YMR142C AND YIL133C AND YNL069C AND YKL006W AND YHL001W AND YLR029C AND YMR121C AND YKL180W AND YJL177W AND YOL120C AND YNL301C AND YMR242C AND YOR312C AND YBR084C-A AND YBL027W AND YBR191W AND YPL079W AND YLR061W AND YFL034C-A AND YBL087C AND YER117W AND YOL127W AND YGL031C AND YGR148C AND YLR344W AND YGR034W AND YHR010W AND YDR471W AND YGL103W AND YFR032C-A AND YGL030W AND YDL075W AND YLR406C AND YBL092W AND YER056C-A AND YIL052C AND YDL191W AND YDL136W AND YPL143W AND YOR234C AND YMR194W AND YPL249C-A AND YNL162W AND YHR141C AND YLR185W AND YDR500C AND YPR043W AND YJR094W-A AND YLR325C AND YJL189W AND YIL148W AND YKR094C AND YDL184C AND YDL133C-A AND YLR340W AND YDL081C AND YDL130W AND YOL039W';
    %% set molecular weight of amino acids [kDa]
    % useful for the weighted sum of p0 and thus for computing proteinWeights
    aaMolW = [89,121,133,147,165,75,155,131,146,131,149,132,115,146,174,105,119,117,204,181];
    aaMolW = aaMolW/1000;
    
    %% formulate enzyme production reactions
    for i=1:length(model.grRules)
        model.grRules{i} = strrep(model.grRules{i},'(','');
        model.grRules{i} = strrep(model.grRules{i},')','');
    end
    
    % parse gene association to figure out how enzymes are built
    [~,~,geneCombs] = parseGeneAssociation(model);

    %% build enzyme production reactions
    newModel = model;    
    newModel.spontaneousRxn = zeros(length(model.rxns),1);
    % compute how many reactions I will end up with and how many enzymes and
    % allocate rxn-enzyme lookup matrix
    nrxnNew = length(model.rxns)+length(unique(geneCombs))+2+1;
    nenzNew = nrxnNew-length(model.rxns);

    newModel.rxnEnzRules = sparse(nrxnNew,nenzNew);
    newModel.enz = {};
	newModel.enzLength = zeros(1,nenzNew);
%     newModel.enzGeneAssoc = cell(nenzNew,1);
    newModel.proteinWeights = zeros(nenzNew,1);
    
    uniqueGeneAssoc = unique(geneCombs);  
    aux = 0; % count how many enzymes I have added up to now
       
    %% add enzyme production reactions
    for i=1:length(uniqueGeneAssoc)
        genesInvolved = strsplit(uniqueGeneAssoc{i},' AND ');
        rxnsCatalysed = false(length(newModel.rxns),1);
        for jj = 1:length(newModel.rxns)
%             rxnsCatalysed(jj) = ~isempty(strfind(newModel.grRules{jj},uniqueGeneAssoc{i}));
            rxnsCatalysed(jj) = strcmp(newModel.grRules{jj},uniqueGeneAssoc{i});
        end
        rxnsCatalysedNames = newModel.rxns(rxnsCatalysed);
        
        % we need a special case for fumarase and glycerol 3-phosphate dehydrogenase, which is the same enzyme
        % catalysing reactions both in cytosol and mitochondrion, so we have
        % to add it twice
        if ismember('r_0451',rxnsCatalysedNames)
            enzName = 'E_r_0451[c11]';
            rxnName = 'synth_E_r_0451';
            aa = sum(AAcounts(ismember(geneNames,genesInvolved),:),1);
            newModel = addMetabolite(newModel,enzName,enzName);
            newModel = addReaction(newModel,{rxnName,['Production of enzyme that catalyzes ',strjoin(rxnsCatalysedNames,' ')]},...
                {%reactants: amino acids
                's_0404[c03]','s_0542[c03]','s_0432[c03]','s_0748[c03]','s_1314[c03]','s_0757[c03]','s_0832[c03]','s_0847[c03]',...
                's_1099[c03]','s_1077[c03]','s_1148[c03]','s_0430[c03]','s_1379[c03]','s_0747[c03]','s_0428[c03]','s_1428[c03]',...
                's_1491[c03]','s_1561[c03]','s_1527[c03]','s_1533[c03]',...% ATP, H2O (one per hydrolysed phosphate), GTP, 10-formyl-THF next
                's_0434[c03]','s_0803[c03]', 's_0785[c03]', 's_0120[c03]',...% products: enzyme, AMP,PPi,GDP,Pi,H+ (one per hydrolysed phosphate), THF, fMet-tRNA-fmet next
                enzName, 's_0423[c03]', 's_0633[c03]', 's_0739[c03]', 's_1322[c03]', 's_0794[c03]', 's_1487[c03]',...%tRNAs
                's_1582[c03]','s_1589[c03]','s_1587[c03]','s_1591[c03]','s_1604[c03]','s_1593[c03]','s_1594[c03]','s_1596[c03]',...
                's_1600[c03]','s_1598[c03]','s_1602[c03]','s_1585[c03]','s_1606[c03]','s_1590[c03]','s_1583[c03]','s_1607[c03]',...
                's_1608[c03]','s_1614[c03]','s_1610[c03]','s_1612[c03]'},...
                [-aa, -sum(aa),-3*sum(aa), -2*sum(aa), -1, 1, sum(aa), sum(aa), 2*sum(aa), 2*sum(aa), 3*sum(aa), 1,aa], false,geneAssocRib);
            newModel.spontaneousRxn = [newModel.spontaneousRxn;0];
            % set rxn-enzyme association
            aux = aux + 1;
            %assert(length(newModel.rxns)==aux+length(model.rxns))
            newModel.rxnEnzRules(findRxnIDs(newModel,'r_0451'),aux) = 1;
            newModel.enz{length(newModel.enz)+1} = enzName;
            newModel.enzLength(aux) = sum(aa);
            %newModel.enzGeneAssoc{length(newModel.enz)} = strjoin(model.genes(genesInvolved));
            newModel.proteinWeights(aux) = sum(aaMolW.*aa);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            enzName = 'E_r_0452[c03]';
            rxnName = 'synth_E_r_0452';
            newModel = addMetabolite(newModel,enzName,enzName);
            newModel = addReaction(newModel,{rxnName,['Production of enzyme that catalyzes ',strjoin(rxnsCatalysedNames,' ')]},...
                {%reactants: amino acids
                's_0404[c03]','s_0542[c03]','s_0432[c03]','s_0748[c03]','s_1314[c03]','s_0757[c03]','s_0832[c03]','s_0847[c03]',...
                's_1099[c03]','s_1077[c03]','s_1148[c03]','s_0430[c03]','s_1379[c03]','s_0747[c03]','s_0428[c03]','s_1428[c03]',...
                's_1491[c03]','s_1561[c03]','s_1527[c03]','s_1533[c03]',...% ATP, H2O,GTP, Met-tRNA-fmet, 10-formyl-THF next
                's_0434[c03]','s_0803[c03]', 's_0785[c03]', 's_0120[c03]',...% products: enzyme, AMP,PPi,GDP,Pi,H+, THF, fMet-tRNA-fmet next
                enzName, 's_0423[c03]', 's_0633[c03]', 's_0739[c03]', 's_1322[c03]', 's_0794[c03]', 's_1487[c03]',...%tRNAs
                's_1582[c03]','s_1589[c03]','s_1587[c03]','s_1591[c03]','s_1604[c03]','s_1593[c03]','s_1594[c03]','s_1596[c03]',...
                's_1600[c03]','s_1598[c03]','s_1602[c03]','s_1585[c03]','s_1606[c03]','s_1590[c03]','s_1583[c03]','s_1607[c03]',...
                's_1608[c03]','s_1614[c03]','s_1610[c03]','s_1612[c03]'},...
                [-aa, -sum(aa), -3*sum(aa),-2*sum(aa), -1, 1, sum(aa), sum(aa), 2*sum(aa), 2*sum(aa),3*sum(aa), 1,aa], false,geneAssocRib);
            newModel.spontaneousRxn = [newModel.spontaneousRxn;0];
            % set rxn-enzyme association
            aux = aux + 1;
            %assert(length(newModel.rxns)==aux+length(model.rxns))
            newModel.rxnEnzRules(findRxnIDs(newModel,'r_0452'),aux) = 1;
            newModel.enz{length(newModel.enz)+1} = enzName;
            newModel.enzLength(aux) = sum(aa);
            %newModel.enzGeneAssoc{length(newModel.enz)} = strjoin(model.genes(genesInvolved));
            newModel.proteinWeights(aux) = sum(aaMolW.*aa);
        elseif ismember('r_0492',rxnsCatalysedNames)
            enzName = 'E_r_0491_1[c03]';
            rxnName = 'synth_E_r_0491_1';
            aa = sum(AAcounts(ismember(geneNames,genesInvolved),:),1);
            newModel = addMetabolite(newModel,enzName,enzName);
            newModel = addReaction(newModel,{rxnName,['Production of enzyme that catalyzes ',strjoin(rxnsCatalysedNames,' ')]},...
                {%reactants: amino acids
                's_0404[c03]','s_0542[c03]','s_0432[c03]','s_0748[c03]','s_1314[c03]','s_0757[c03]','s_0832[c03]','s_0847[c03]',...
                's_1099[c03]','s_1077[c03]','s_1148[c03]','s_0430[c03]','s_1379[c03]','s_0747[c03]','s_0428[c03]','s_1428[c03]',...
                's_1491[c03]','s_1561[c03]','s_1527[c03]','s_1533[c03]',...% ATP, H2O (one per hydrolysed phosphate), GTP, 10-formyl-THF next
                's_0434[c03]','s_0803[c03]', 's_0785[c03]', 's_0120[c03]',...% products: enzyme, AMP,PPi,GDP,Pi,H+ (one per hydrolysed phosphate), THF, fMet-tRNA-fmet next
                enzName, 's_0423[c03]', 's_0633[c03]', 's_0739[c03]', 's_1322[c03]', 's_0794[c03]', 's_1487[c03]',...%tRNAs
                's_1582[c03]','s_1589[c03]','s_1587[c03]','s_1591[c03]','s_1604[c03]','s_1593[c03]','s_1594[c03]','s_1596[c03]',...
                's_1600[c03]','s_1598[c03]','s_1602[c03]','s_1585[c03]','s_1606[c03]','s_1590[c03]','s_1583[c03]','s_1607[c03]',...
                's_1608[c03]','s_1614[c03]','s_1610[c03]','s_1612[c03]'},...
                [-aa, -sum(aa),-3*sum(aa), -2*sum(aa), -1, 1, sum(aa), sum(aa), 2*sum(aa), 2*sum(aa), 3*sum(aa), 1,aa], false,geneAssocRib);
            newModel.spontaneousRxn = [newModel.spontaneousRxn;0];
            % set rxn-enzyme association
            aux = aux + 1;
            %assert(length(newModel.rxns)==aux+length(model.rxns))
            newModel.rxnEnzRules(findRxnIDs(newModel,'r_0491_1'),aux) = 1;
            newModel.enz{length(newModel.enz)+1} = enzName;
            newModel.enzLength(aux) = sum(aa);
            %newModel.enzGeneAssoc{length(newModel.enz)} = strjoin(model.genes(genesInvolved));
            newModel.proteinWeights(aux) = sum(aaMolW.*aa);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            enzName = 'E_r_0492[c11]';
            rxnName = 'synth_E_r_0492';
            newModel = addMetabolite(newModel,enzName,enzName);
            newModel = addReaction(newModel,{rxnName,['Production of enzyme that catalyzes ',strjoin(rxnsCatalysedNames,' ')]},...
                {%reactants: amino acids
                's_0404[c03]','s_0542[c03]','s_0432[c03]','s_0748[c03]','s_1314[c03]','s_0757[c03]','s_0832[c03]','s_0847[c03]',...
                's_1099[c03]','s_1077[c03]','s_1148[c03]','s_0430[c03]','s_1379[c03]','s_0747[c03]','s_0428[c03]','s_1428[c03]',...
                's_1491[c03]','s_1561[c03]','s_1527[c03]','s_1533[c03]',...% ATP, H2O,GTP, Met-tRNA-fmet, 10-formyl-THF next
                's_0434[c03]','s_0803[c03]', 's_0785[c03]', 's_0120[c03]',...% products: enzyme, AMP,PPi,GDP,Pi,H+, THF, fMet-tRNA-fmet next
                enzName, 's_0423[c03]', 's_0633[c03]', 's_0739[c03]', 's_1322[c03]', 's_0794[c03]', 's_1487[c03]',...%tRNAs
                's_1582[c03]','s_1589[c03]','s_1587[c03]','s_1591[c03]','s_1604[c03]','s_1593[c03]','s_1594[c03]','s_1596[c03]',...
                's_1600[c03]','s_1598[c03]','s_1602[c03]','s_1585[c03]','s_1606[c03]','s_1590[c03]','s_1583[c03]','s_1607[c03]',...
                's_1608[c03]','s_1614[c03]','s_1610[c03]','s_1612[c03]'},...
                [-aa, -sum(aa), -3*sum(aa),-2*sum(aa), -1, 1, sum(aa), sum(aa), 2*sum(aa), 2*sum(aa),3*sum(aa), 1,aa], false,geneAssocRib);
            newModel.spontaneousRxn = [newModel.spontaneousRxn;0];
            % set rxn-enzyme association
            aux = aux + 1;
            %assert(length(newModel.rxns)==aux+length(model.rxns))
            newModel.rxnEnzRules(findRxnIDs(newModel,'r_0492'),aux) = 1;
            newModel.enz{length(newModel.enz)+1} = enzName;
            newModel.enzLength(aux) = sum(aa);
            %newModel.enzGeneAssoc{length(newModel.enz)} = strjoin(model.genes(genesInvolved));
            newModel.proteinWeights(aux) = sum(aaMolW.*aa);
        else
            [~, tmp_met_struct] = regexp(newModel.mets,'(?<met>.+)\[(?<comp>.+)\]|(?<met>.+)\((?<comp>.+)\)','tokens','names');
            tmp_met_struct = [tmp_met_struct{:}];
            ridx = findRxnIDs(newModel,rxnsCatalysedNames);
            midx = [];
            for l=1:length(rxnsCatalysedNames)
                midx = [midx;find(model.S(:,ridx(l)))];
            end
            comps = unique({tmp_met_struct(midx).comp});
            ending = '';
            for j=1:length(comps)
                ending = [ending,'_[',comps{j},']'];
            end
            ending = ending(2:end);
    
            enzName = sprintf('E_%s%s',strjoin(rxnsCatalysedNames,'_'),ending);
            if any(strncmp(newModel.enz,enzName,length(enzName)))
                enzName = sprintf('%s_%d%s',enzName,sum(strncmp(newModel.enz,enzName,length(enzName)))+1,ending);
            end
            
            % make sure we take into account possible stoichiometry
            % information for the enzymes
            idx = find(ismember(T.geneAssoc,uniqueGeneAssoc{i}));
            if ~isempty(idx)
                stoich = str2num(table2array(T.stoichiometry(idx)));
                aa = zeros(1,20);
                for d=1:length(stoich)
                    aa = aa + stoich(d)*AAcounts(ismember(geneNames,genesInvolved(d)),:);
                end
            end
            
            newModel = addMetabolite(newModel,enzName,enzName);
            rxnName = strcat('synth_',sprintf('E_%s%s',strjoin(rxnsCatalysedNames,'_')));
            newModel = addReaction(newModel,{rxnName,['Production of enzyme that catalyzes ',strjoin(rxnsCatalysedNames,' ')]},...
                {%reactants: amino acids
                's_0404[c03]','s_0542[c03]','s_0432[c03]','s_0748[c03]','s_1314[c03]','s_0757[c03]','s_0832[c03]','s_0847[c03]',...
                's_1099[c03]','s_1077[c03]','s_1148[c03]','s_0430[c03]','s_1379[c03]','s_0747[c03]','s_0428[c03]','s_1428[c03]',...
                's_1491[c03]','s_1561[c03]','s_1527[c03]','s_1533[c03]',...% ATP, GTP, 10-formyl-THF next
                's_0434[c03]','s_0803[c03]', 's_0785[c03]', 's_0120[c03]',...% products: enzyme, AMP,PPi,GDP,Pi, THF, fMet-tRNA-fmet next
                enzName, 's_0423[c03]', 's_0633[c03]', 's_0739[c03]', 's_1322[c03]','s_0794[c03]', 's_1487[c03]',...%tRNAs
                's_1582[c03]','s_1589[c03]','s_1587[c03]','s_1591[c03]','s_1604[c03]','s_1593[c03]','s_1594[c03]','s_1596[c03]',...
                's_1600[c03]','s_1598[c03]','s_1602[c03]','s_1585[c03]','s_1606[c03]','s_1590[c03]','s_1583[c03]','s_1607[c03]',...
                's_1608[c03]','s_1614[c03]','s_1610[c03]','s_1612[c03]'},...
                [-aa, -sum(aa),-3*sum(aa), -2*sum(aa), -1, 1, sum(aa), sum(aa), 2*sum(aa), 2*sum(aa),3*sum(aa), 1,aa], false, geneAssocRib);
            newModel.spontaneousRxn = [newModel.spontaneousRxn;0];
            % set rxn-enzyme association
            aux = aux + 1;
            newModel.rxnEnzRules(rxnsCatalysed,aux) = 1;
            newModel.enz{length(newModel.enz)+1} = enzName;
            newModel.enzLength(aux) = sum(aa);      
            newModel.proteinWeights(aux) = sum(aaMolW.*aa);
        end        
    end
    % check everything went well
    assert(length(newModel.enz)==aux)
    assert(nrxnNew==length(model.rxns)+aux+1)
 
    %% mark spontaneous reactions
    for i=1:length(model.grRules)
        if strcmp(model.grRules{i},'')
            newModel.spontaneousRxn(i) = 1;
        end
    end

    %% make ribosome
    % gene composition from KEGG + ribosomal RNA
    M = importdata('AminoAcidCountsRibosomeYeast.csv',',',1);
    aa = M.data;
    M = importdata('ribosomalRNAcounts.csv',',',1);
    BaseCountsRNA = M.data(end,:);
    % add energy requirement for production of ribosome into BaseCountsRNA
    BaseCountsRNA(1) = BaseCountsRNA(1) + sum(aa);
	BaseCountsRNA(3) = BaseCountsRNA(3) + 2*sum(aa);
    enzName = 'Ribosome[c03]';
    newModel = addMetabolite(newModel,enzName,'Ribosome[c03]');
    newModel = addReaction(newModel,{'synth_Ribosome','Production of ribosome'},...
        {%reactants: amino acids
        's_0404[c03]','s_0542[c03]','s_0432[c03]','s_0748[c03]','s_1314[c03]','s_0757[c03]','s_0832[c03]','s_0847[c03]',...
        's_1099[c03]','s_1077[c03]','s_1148[c03]','s_0430[c03]','s_1379[c03]','s_0747[c03]','s_0428[c03]','s_1428[c03]',...
        's_1491[c03]','s_1561[c03]','s_1527[c03]','s_1533[c03]',...% ATP, CTP, GTP, UTP, 10-formyl-THF,h2O next
        's_0434[c03]','s_0539[c03]','s_0785[c03]','s_1559[c03]','s_0120[c03]','s_0803[c03]',...% products: enzyme, AMP,PPi,GDP,Pi, THF, fMet-tRNA-fmet, h+ next
        'Ribosome[c03]','s_0423[c03]', 's_0633[c03]', 's_0739[c03]', 's_1322[c03]', 's_1487[c03]','s_0794[c03]',...%tRNAs
        's_1582[c03]','s_1589[c03]','s_1587[c03]','s_1591[c03]','s_1604[c03]','s_1593[c03]','s_1594[c03]','s_1596[c03]',...
        's_1600[c03]','s_1598[c03]','s_1602[c03]','s_1585[c03]','s_1606[c03]','s_1590[c03]','s_1583[c03]','s_1607[c03]',...
        's_1608[c03]','s_1614[c03]','s_1610[c03]','s_1612[c03]'},...
        [-aa,-BaseCountsRNA,-1,-3*sum(aa),1,sum(aa),sum(aa)+sum(BaseCountsRNA),2*sum(aa),2*sum(aa),1,3*sum(aa),aa],false,geneAssocRib);
    newModel.spontaneousRxn = [newModel.spontaneousRxn;0];
    % ribosome catalyses all protein formation reactions including its own
    newModel.rxnEnzRules(length(model.rxns)+1:end,end) = 1;
    
    newModel.proteinWeights(end) = sum(aaMolW.*aa);
    newModel.enz{length(newModel.enz)+1} = 'Ribosome[c03]';
    newModel.sizePmet = newModel.sizeQuotaMet + length(newModel.enz);
    newModel.sizePrxn = newModel.sizeQuotaRxn + nrxnNew-length(model.rxns);   
	newModel.enzLength(end) = sum(aa);
    %newModel.enzGeneAssoc{length(newModel.enz)} = 'Ribosome';
end
