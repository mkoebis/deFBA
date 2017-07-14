function sbmlModel = writeSBML_ram(model,fileName,compSymbolList,compNameList, storageCompID, biomassCompID,externalCompID)

% writeSBML_ram exports a deFBA model structure into an SBML FBCv2 file using the ram standard
%
%INPUTS
% model             deFBA model structure with the fields:
%   rxns                        cell array of all reaction IDs
%   rxnNames                    cell array of all metabolite names
%   rev                         0-1 array describing if the reactions are reversible (1) or irreversible (0)
%   spontaneousRxn              0-1 array describing if the reactions are spontaneous (1) or enzyme-catalyzed (0)
%   noRxn                       number or all reactions
%   grRules                     cell array with strings describing the gene association for each reaction
%   rxnECNumbers                cell array containing the EC number for each reaction
%   mets                        cell array of all species IDs
%   metNames                    cell array of all species names
%   S                           stoichiometric matrix
%   genes                       cell array of all gene names
%   rxnGeneMat                  length(rxns)xlength(genes) 0-1 matrix describing if gene j is involved in the catalysis of reaction i (entry at row i, column j is 1) or not (entry at row i, column j is 0)
%   tf                          final time of simulation (hours)
%   N                           number of discretization points for deFBA timecourses
%   epsilon                     scaling factor for improved numerics (see Waldherr et al 2015)
%   beta                        dilution factor for enzymes, quota and storage
%   sizeXmet                    number of internal metabolites
%   sizeYmet                    number of external and storage metabolites
%   sizeQuotaMet                number of quota components
%   sizePmet                    number of quota and enzyme species
%   sizeXrxn                    number of reactions producing internal metabolites
%   sizeYrxn                    number of reactions involving storage or external metabolites
%   sizeQuotaRxn                number of reactions producing quota metabolites
%   sizePrxn                    number of reactions producing quota or enzymes
%   noStorage                   number of storage metabolites
%   storageWeight               array with molecular weight of storage metabolites in kDa
%   enz                         cell array with species IDs of enzymes
%   proteinWeights              array storing the molecular weights of the enzymes in kDa
%   rxnEnzRules                 length(rxns)xlength(enz) 0-1 matrix describing if enzyme j catalysises reaction i (entry at row i, column j is 1) or not (entry at row i, column j is 0)
%   Kcat_f                      length(rxns)xlength(enz) matrix storing the kcat of enzyme j for the forward direction of reaction i (entry at row i, column j is kcat)
%   quotaInitial                array containing the total biomass percentages that have to be satisfied for each quota compound at each time point
%   initialBiomass              initial values for the quota and enzymes amounts
%   Y0                          initial values for the storage (first noStorage entries) and external metabolites (rest) amounts
%   maintenanceID               ID of the maintenance reaction
%   maintenanceValue            maintenance dependence on total biomass
%   gprComp                     cell array that specifies for each reaction the recipe for building its corresponding enzyme from the individual gene ids ( e.g. 3*gene1 AND 2*gene2 means that the enzyme is made of three copies of gene1 and 2 copies of gene2)
% fileName          File name for output file
% compSymbolList    List of compartment ids
% compNameList      List of compartment names corresponding to compSymbolList
% externalCompID    String containing the ID of the extracellular metabolites compartment in the SBML file
% biomassCompID     String containing the ID of the biomass components compartment in the SBML file
% storageCompID     String containing the ID of the storage components compartment in the SBML file
%
%OUTPUT
% sbmlModel         SBML MATLAB structure
%
% Alexandra Reimers 13/07/2017

%% define empty structs for sbml
sbmlModel = initSBML(model);

%% build temporary structs
tmp_unitDefinition=struct('typecode','SBML_UNIT_DEFINITION',...
        'metaid','',...
        'notes','',...
        'annotation','',...
        'sboTerm',-1,...
        'name','',...
        'id','',...
        'unit','',...
        'level',3,...
        'version',1);
    tmp_unitDefinition.unit =struct('typecode','SBML_UNIT',...
        'metaid','',...
        'notes','',...
        'annotation','',...
        'sboTerm',-1,...
        'kind','',...
        'exponent','',...
        'scale','',...
        'multiplier','',...
        'level',3,...
        'version',1);
    tmp_species=struct('typecode','SBML_SPECIES',...
        'metaid','',...
        'notes','',...
        'annotation','',...
        'sboTerm',-1,...
        'name','',...
        'id','',...  
        'compartment','',... 
        'initialAmount',0,...  
        'initialConcentration','',...
        'substanceUnits','',...
        'hasOnlySubstanceUnits',0,...
        'boundaryCondition',0,...
        'constant',0,...
        'conversionFactor','',...
        'isSetInitialAmount',1,...
        'isSetInitialConcentration','',...
        'fbc_charge','',...             
        'fbc_chemicalFormula','',...  
        'isSetfbc_charge','',...  
        'level',3,...    
        'version',1,...
        'fbc_version',2);
    tmp_parameter=struct('typecode','SBML_PARAMETER',...
        'metaid','',...
        'notes','',...
        'annotation','',...
        'sboTerm',-1,...
        'name','',...
        'id','noID',... 
        'value','',...
        'units','',...
        'constant',1,...  
        'isSetValue',1,...
        'level',3,...
        'version',1);
    tmp_compartment=struct('typecode','SBML_COMPARTMENT',...
        'metaid','',...
        'notes','',...
        'annotation','',...
        'sboTerm',-1,...
        'name','',... 
        'id','',...     
        'spatialDimensions',3,... 
        'size',1,...
        'units','',...
        'constant',1,...
        'isSetSize',1,...
        'isSetSpatialDimensions',1,...
        'level',3,...
        'version',1);
    tmp_Rxn=struct('typecode','SBML_REACTION',...
        'metaid','',...
        'notes','',... 
        'annotation','',...
        'sboTerm',-1,...
        'name','',... 
        'id','',... 
        'reactant','',...
        'product','',...
        'modifier','',...
        'kineticLaw','',...
        'reversible','',... 
        'fast','',...
        'compartment','',...
        'isSetFast','',...
        'fbc_lowerFluxBound','',...
        'fbc_upperFluxBound','',...
        'fbc_geneProductAssociation','',...
        'level',3,...
        'version',1,...
        'fbc_version',2);
    sbml_tmp_species_ref=struct('typecode','SBML_SPECIES_REFERENCE',... 
        'metaid','',...
        'notes','',...
        'annotation','',...
        'sboTerm',-1,...
        'species','',... 
        'id','',...
        'name','',...
        'stoichiometry','',... 
        'constant','',...
        'isSetStoichiometry','',...
        'level',3,...
        'version',1);
    tmp_Rxn.modifier=struct('typecode',cell(1,0),...
        'metaid',cell(1,0),...
        'notes',cell(1,0),...
        'annotation',cell(1,0),...
        'sboTerm',-1,...
        'species',cell(1,0),...
        'id',cell(1,0),...
        'name',cell(1,0),...
        'level',3,...
        'version',1);
    tmp_Rxn.kineticLaw=struct('typecode',cell(1,0),...
        'metaid',cell(1,0),...
        'notes',cell(1,0),...
        'annotation',cell(1,0),...
        'sboTerm',-1,...
        'math',cell(1,0),...
        'localParameter',cell(1,0),...
        'level',3,...
        'version',1);
    tmp_fbc_fluxBound=struct('typecode','SBML_FBC_FLUXBOUND',...
        'metaid','',...
        'notes','',...
        'annotation','',...
        'sboTerm',-1,...
        'fbc_id','',...
        'fbc_reaction','',... 
        'fbc_operation','',...
        'fbc_value','',...    
        'isSetfbc_value','',...
        'level',3,...
        'version',1,...
        'fbc_version',2);
    tmp_Rxn.fbc_geneProductAssociation=struct('typecode','SBML_FBC_GENE_PRODUCT_ASSOCIATION',...
        'metaid','',...
        'notes','',...
        'annotation','',...
        'sboTerm',-1,...
        'fbc_id','',...
        'fbc_name','',...
        'fbc_association','',...
        'level',3,...
        'version',1,...
        'fbc_version',2);
    tmp_Rxn.fbc_geneProductAssociation.fbc_association=struct('typecode','SBML_FBC_GENE_PRODUCT_REF',...
        'metaid','',... 
        'notes','',...
        'annotation','',...
        'sboTerm',-1,...
        'fbc_association','',...
        'level',3,...
        'version',1,...
        'fbc_version',2);
    tmp_fbc_geneProduct=struct('typecode','SBML_FBC_GENE_PRODUCT',...
        'metaid','',... 
        'notes','',...
        'annotation','',...
        'sboTerm',-1,...
        'fbc_id','',...
        'fbc_name','',...
        'fbc_label','',...
        'fbc_associatedSpecies','',...
        'level',3,...
        'version',1,...
        'fbc_version',2);

%% units
reaction_units = 'mmol_per_hr';
tmp_unitDefinition.id =  reaction_units;
unit_kinds = {'mole','gram','second'};
unit_exponents = [1 -1 -1];
unit_scales = [-3 0 0];
unit_multipliers = [1 1 1*60*60];
sbmlModel.unitDefinition=tmp_unitDefinition;

for i = 1:size(unit_kinds, 2)
    tmp_unitDefinition.unit.kind = unit_kinds{ i };
    tmp_unitDefinition.unit.exponent = unit_exponents(i);
    tmp_unitDefinition.unit.scale = unit_scales(i);
    tmp_unitDefinition.unit.multiplier = unit_multipliers(i);
    if i==1
        sbmlModel.unitDefinition.unit=tmp_unitDefinition.unit;
    else
        sbmlModel.unitDefinition.unit=[sbmlModel.unitDefinition.unit,tmp_unitDefinition.unit];
    end
end

%% species
tmp_metCompartment = {};

[~, tmp_met_struct] = regexp(model.mets,'(?<met>.+)\[(?<comp>.+)\]|(?<met>.+)\((?<comp>.+)\)','tokens','names');

sbmlModel.species=tmp_species;

% already create parameter for zero
tmp_parameter.id = 'zero';
tmp_parameter.value = 0;
sbmlModel.parameter = tmp_parameter;

%% metabolites
for i=1:length(model.mets)
    tmp_note='';
    if  ~isempty(tmp_met_struct{1})
        tmp_met = tmp_met_struct{i}.met;
    else
        tmp_met = model.mets{i};
    end
    
    if isempty(tmp_met_struct{1})
        %Change id to correspond to SBML id specifications
        tmp_met = strcat((tmp_met), '_', '[c]');
    else
        tmp_met = strcat((tmp_met), '_', tmp_met_struct{i}.comp);
    end
    
    tmp_species.id = formatID(tmp_met);
    tmp_species.name = model.metNames{i};
    tmp_species.compartment = formatID(tmp_met_struct{i}.comp);
    tmp_metCompartment{i} = formatID(tmp_met_struct{i}.comp); 


%% annotations for metabolites
 tmp_noteBegin = sprintf('<annotation>\n<ram:RAM xmlns:ram="http://fancy.web.page">\n<ram:species');
    if strcmp(tmp_metCompartment{i},biomassCompID)        
        if ~isempty(strfind(tmp_species.name,'quota'))
            [~,idxQ] = ismember(tmp_species.name,model.metNames(model.sizeXmet+model.sizeYmet+1:model.sizeXmet+model.sizeYmet+model.sizeQuotaMet));
            tmp_percBiomass = model.quotaInitial(idxQ);
            tmp_parameter.id = sprintf('bioPercentage_%d',idxQ);
            idP = tmp_parameter.id;
            tmp_parameter.value = tmp_percBiomass;
            sbmlModel.parameter = [sbmlModel.parameter,tmp_parameter];
            tmp_molWeight = model.quotaWeights(idxQ);
            tmp_parameter.id = sprintf('weight_quota_%d',idxQ);
            idW = tmp_parameter.id;
            tmp_parameter.value = tmp_molWeight;
            sbmlModel.parameter = [sbmlModel.parameter,tmp_parameter];
            tmp_note=sprintf('%s ram:molecularWeight="%s" ram:biomassPercentage="%s"/>\n</ram:RAM>\n</annotation>',tmp_noteBegin,idW,idP);
        else
            [~,idxP] = ismember(tmp_species.name,model.metNames(model.sizeXmet+model.sizeYmet+model.sizeQuotaMet+1:model.sizeXmet+model.sizeYmet+model.sizePmet));
            tmp_molWeight = model.proteinWeights(idxP);
            tmp_parameter.id = sprintf('weight_%d',idxP);
            tmp_parameter.value = tmp_molWeight;
            sbmlModel.parameter = [sbmlModel.parameter,tmp_parameter];
            tmp_note=sprintf('%s ram:molecularWeight="%s" ram:biomassPercentage="zero"/>\n</ram:RAM>\n</annotation>',tmp_noteBegin,tmp_parameter.id);
        end
    end

    if strcmp(tmp_metCompartment{i},storageCompID)
        [~,idxY] = ismember(tmp_species.name,model.metNames(model.sizeXmet+1:model.sizeXmet+model.sizeYmet));
        tmp_molWeight = model.storageWeight(idxY);
        tmp_parameter.id = sprintf('weight_st_%d',idxY);
        tmp_parameter.value = tmp_molWeight;
        sbmlModel.parameter = [sbmlModel.parameter,tmp_parameter];
        tmp_note=sprintf('%s ram:molecularWeight="%s" ram:biomassPercentage="zero"/>\n</ram:RAM>\n</annotation>',tmp_noteBegin,tmp_parameter.id);
    end
    
    tmp_species.annotation = tmp_note;
    
    if i==1
        sbmlModel.species=tmp_species;
    else
        sbmlModel.species=[sbmlModel.species, tmp_species];
    end
    
end

%% compartments.

tmp_metCompartment = unique(tmp_metCompartment);

for i=1:size(tmp_metCompartment,2)
    if ~isempty(tmp_metCompartment) % in the case of an empty model
        tmp_id = tmp_metCompartment{1,i};
        tmp_name = compNameList{find(strcmp(formatID(compSymbolList),tmp_id))};  
        
        tmp_compartment.id = formatID(tmp_id);
        tmp_compartment.name = tmp_name;
    end
    
    if i==1
        sbmlModel.compartment=tmp_compartment;
    else
        sbmlModel.compartment=[sbmlModel.compartment, tmp_compartment];
    end
end

%% reactions
sbmlModel.reaction = tmp_Rxn;
sbmlModel.fbc_geneProduct = tmp_fbc_geneProduct; 

%% add gene products with correct stoichiometry
geneProductRxn = cell(length(model.rxns),1);
if isfield(model,'gprComp')
    t = 1;
    for i=1:length(model.gprComp)
        % find index of corresponding enzyme in the model
        idx = find(model.rxnEnzRules(i,:));
        % find index of corresponding enzyme in sbml
        idx = model.sizeXmet+model.sizeYmet+model.sizeQuotaMet+idx;        
        idx2 = find(model.rxnEnzRules(i,:));
             
        % correct ids and names of enzymes
        gIdx = find(model.rxnGeneMat(i,:));
        if length(gIdx)==1            
            rIdx = find(model.rxnEnzRules(:,idx2));
            for j=1:length(rIdx)
                if j==1
                    aux = sprintf('enzyme catalyzing reaction(s) %s',model.rxns{rIdx(j)});
                else
                    aux = sprintf('%s %s',aux,model.rxns{rIdx(j)});
                end
            end     
            aux = sprintf('%s[bio]',aux);
            sbmlModel.species(idx).name = aux;
            % find metabolites involved
            [mIdx,~] = find(model.S(:,rIdx));
            %get their compartments
            comps = unique({sbmlModel.species(mIdx).compartment});
            sbmlModel.species(idx).id = model.genes{gIdx};
            % add the compartments to the name
            for j=1:length(comps)
                sbmlModel.species(idx).id = [sbmlModel.species(idx).id,'_',comps{j}];
            end
            sbmlModel.species(idx).id = [sbmlModel.species(idx).id,'_bio'];            
            geneProductRxn{i} = sbmlModel.species(idx).id;
        elseif length(gIdx)>1
            if isempty(strfind(sbmlModel.species(idx).id,'complex'))
                rIdx = find(model.rxnEnzRules(:,idx2));
                for j=1:length(rIdx)
                    if j==1
                        aux = sprintf('enzyme complex catalyzing reaction(s) %s',model.rxns{rIdx(j)});
                    else
                        aux = sprintf('%s %s',aux,model.rxns{rIdx(j)});
                    end
                end     
                aux = sprintf('%s[bio]',aux);
                sbmlModel.species(idx).name = aux;
                
                 % find metabolites involved
                [mIdx,~] = find(model.S(:,rIdx));
                %get their compartments
                comps = unique({sbmlModel.species(mIdx).compartment});
                sbmlModel.species(idx).id = sprintf('complex_%d',t);
                t = t+1;
                % add the compartments to the name
                for j=1:length(comps)
                    sbmlModel.species(idx).id = [sbmlModel.species(idx).id,'_',comps{j}];
                end
                sbmlModel.species(idx).id = [sbmlModel.species(idx).id,'_bio'];
                
                geneProductRxn{i} = sbmlModel.species(idx).id;
            else
                geneProductRxn{i} = sbmlModel.species(idx).id;
            end
        end
        
    end
end

% set initial amounts for metabolites allowed to accumulate and biomass
t = 1;
b = 1;
c = 1;
for i=1:length(sbmlModel.species)
   
    if strcmp(sbmlModel.species(i).compartment,storageCompID)
        sbmlModel.species(i).initialAmount = model.Y0(t);
        sbmlModel.species(i).constant = 0;
        t = t+1;
    elseif strcmp(sbmlModel.species(i).compartment,externalCompID)
        if ~ismember(sbmlModel.species(i).id,model.carbonSources)
            sbmlModel.species(i).initialAmount = 0;
            sbmlModel.species(i).constant = 1;
            sbmlModel.species(i).boundaryCondition = 1;
        else
            sbmlModel.species(i).initialAmount = model.Y0(length(model.storageWeight)+c);
            sbmlModel.species(i).constant = 0;
            sbmlModel.species(i).boundaryCondition = 1;
            c = c+1;
        end
    elseif strcmp(sbmlModel.species(i).compartment,biomassCompID)
        sbmlModel.species(i).initialAmount = model.initialBiomass(b);
        b = b+1;
    end
end

gprIdx = true(length(geneProductRxn),1);
for i=1:length(geneProductRxn)
    if isempty(geneProductRxn{i})
        gprIdx(i) = false;
    end
end
geneProductRxnOrig = geneProductRxn;
gprComp = model.gprComp(gprIdx);
geneProductRxn = geneProductRxn(gprIdx);
[geneProductRxn,idxUnique,~] = unique(geneProductRxn);

gprComp = gprComp(idxUnique);
[~,u] = unique(gprComp);
doubles = gprComp(setdiff(1:1:length(gprComp),u));
countDummy=1;
for i=1:length(geneProductRxn)
    tmp_fbc_geneProduct.fbc_id=formatID(geneProductRxn{i}); % generate fbc_id values
    if ismember(gprComp{i},doubles)
        tmp_fbc_geneProduct.fbc_label=sprintf('%s_%d',gprComp{i},countDummy);
        countDummy = countDummy+1;
    else
        tmp_fbc_geneProduct.fbc_label=gprComp{i};
    end
    tmp_fbc_geneProduct.fbc_associatedSpecies = formatID(geneProductRxn{i});
    if i==1
        sbmlModel.fbc_geneProduct=tmp_fbc_geneProduct;
    else
        sbmlModel.fbc_geneProduct=[sbmlModel.fbc_geneProduct,tmp_fbc_geneProduct];
    end
end

% add maintenance scaling parameter
tmp_parameter.id = 'maintenancePercentage';
tmp_parameter.value = model.maintenanceValue;
sbmlModel.parameter = [sbmlModel.parameter,tmp_parameter];
count=1;
for i=1:size(model.rxns, 1)   
    % write EC number
    t = tmp_Rxn;   
    
    if isfield(model,'rxnECNumbers')&&i<=length(model.rxnECNumbers)
        tmp_note = [' <p>EC Number: ' model.rxnECNumbers{i} '</p>'];
    end
    tmp_note = ['<body xmlns="http://www.w3.org/1999/xhtml">' tmp_note '</body>'];
    t.notes=tmp_note;
      
    tmp_rxnID =  formatID(model.rxns{i});
    
    tmp_rxnName='';
    if isfield(model, 'rxnNames')
        tmp_rxnName = model.rxnNames{i};
    end
    
    tmp_rxnRev= model.rev(i);    
    kcat = model.Kcat_f(i,(model.Kcat_f(i,:)~=0));
    isMaintenance = strcmp(model.rxns{i},model.maintenanceID);
    tmp_note = sprintf('<annotation>\n<ram:RAM xmlns:ram="http://fancy.web.page">\n<ram:reaction');
    if isMaintenance
        tmp_note = sprintf('%s ram:maintenanceScaling="maintenancePercentage" ram:kcatForward="zero" ram:kcatBackward="zero"/>\n</ram:RAM>\n</annotation>',tmp_note);
    else
        % add kcat to list of parameters and then link
        if model.spontaneousRxn(i)==0
            tmp_parameter.id = sprintf('kcat_%d',count);
            tmp_parameter.value = kcat;
            count = count + 1;
            sbmlModel.parameter = [sbmlModel.parameter,tmp_parameter];
            if model.rev(i)
                tmp_note = sprintf('%s ram:maintenanceScaling="zero" ram:kcatForward="%s" ram:kcatBackward="%s"/>\n</ram:RAM>\n</annotation>',tmp_note,tmp_parameter.id,tmp_parameter.id);
            else
                tmp_note = sprintf('%s ram:maintenanceScaling="zero" ram:kcatForward="%s" ram:kcatBackward="zero"/>\n</ram:RAM>\n</annotation>',tmp_note,tmp_parameter.id);
            end
        else
            tmp_note = sprintf('%s ram:maintenanceScaling="zero" ram:kcatForward="zero" ram:kcatBackward="zero"/>\n</ram:RAM>\n</annotation>',tmp_note);
        end
    end
    t.annotation = tmp_note;
        
    if i <= model.noRxn-model.sizePrxn
        t.fast = 1;
    else
        t.fast = 0;
    end
    t.id=tmp_rxnID;
    t.name=tmp_rxnName;
    t.reversible=tmp_rxnRev;    
    t.isSetFast=1;
    
    %Add in the reactants and products
    met_idx = find(model.S(:, i));
    t.product=[];
    t.reactant=[];
    for j_met=1:size(met_idx,1)
        tmp_idx = met_idx(j_met,1);
        sbml_tmp_species_ref.species = sbmlModel.species(tmp_idx).id; % model.mets{tmp_idx};
        met_stoich = model.S(tmp_idx, i);
        sbml_tmp_species_ref.stoichiometry = abs(met_stoich);
        sbml_tmp_species_ref.isSetStoichiometry=1;
        sbml_tmp_species_ref.constant=1;
        if (met_stoich > 0)
            t.product = [ t.product, sbml_tmp_species_ref ];
        else
            t.reactant = [ t.reactant, sbml_tmp_species_ref];
        end
    end
    %% grRules
        if ~isempty(geneProductRxnOrig{i})
            t.fbc_geneProductAssociation.fbc_id=sprintf('ga_%d',i);
            t.fbc_geneProductAssociation.fbc_association.fbc_association=formatID(geneProductRxnOrig{i});
        end
   
    if i==1;
        sbmlModel.reaction=t;
    else
        sbmlModel.reaction=[sbmlModel.reaction,t];
    end
    
end

sbmlModel.fbc_strict=0; 
sbmlModel.fbc_objective='';

fbcStr=['http://www.sbml.org/sbml/level3/version1/','fbc/version',num2str(2)];

sbmlModel.namespaces=struct('prefix',{'','fbc'},...
    'uri',{'http://www.sbml.org/sbml/level3/version1/core',...
    fbcStr});

OutputSBML(sbmlModel,fileName);
end

%% initialize SBML model
function sbmlModel = initSBML(model)
    sbmlModel.constraint=struct('typecode',cell(1,0),...
        'metaid',cell(1,0),...
        'notes',cell(1,0),...
        'annotation','',...
        'sboTerm',-1,...
        'math',cell(1,0),...
        'message',cell(1,0),...
        'level',3,...
        'version',1);
    sbmlModel.functionDefinition=struct('typecode',cell(1,0),...
        'metaid',cell(1,0),...
        'notes',cell(1,0),...
        'annotation',cell(1,0),...
        'sboTerm',-1,...
        'name',cell(1,0),...
        'id',cell(1,0),...
        'math',cell(1,0),...
        'level',3,...
        'version',1);
    sbmlModel.event=struct('typecode',cell(1,0),...
        'metaid',cell(1,0),...
        'notes',cell(1,0),...
        'annotation',cell(1,0),...
        'sboTerm',-1,...
        'name',cell(1,0),...
        'id',cell(1,0),...
        'useValuesFromTriggerTime',cell(1,0),...
        'trigger',cell(1,0),...
        'delay',cell(1,0),...
        'priority',cell(1,0),...
        'eventAssignment',cell(1,0),...
        'level',3,...
        'version',1);
    sbmlModel.rule=struct('typecode',cell(1,0),...
        'metaid',cell(1,0),...
        'notes',cell(1,0),...
        'annotation',cell(1,0),...
        'sboTerm',-1,...
        'formula',cell(1,0),...
        'variable',cell(1,0),...
        'species',cell(1,0),...
        'compartment',cell(1,0),...
        'name',cell(1,0),...
        'units',cell(1,0),...
        'level',3,...
        'version',1);
    sbmlModel.initialAssignment=struct('typecode',cell(1,0),...
        'metaid',cell(1,0),...
        'notes',cell(1,0),...
        'annotation',cell(1,0),...
        'sboTerm',-1,...
        'symbol',cell(1,0),...
        'math',cell(1,0),...
        'level',3,...
        'version',1);
    list={'SBML_level';
        'SBML_version';
        'annotation';
        'areaUnits';
        'avogadro_symbol';
        'conversionFactor';
        'delay_symbol';
        'extentUnits';
        'fbc_activeObjective';
        'fbc_version';
        'id';
        'lengthUnits';
        'metaid';
        'name';
        'notes';
        'sboTerm';
        'substanceUnits';
        'timeUnits';
        'time_symbol';
        'typecode';
        'volumeUnits'};
    for i=1:length(list)
        if strfind(list{i},'SBML_level')
            sbmlModel.(list{i})=3;
        elseif strfind(list{i},'SBML_version')
            sbmlModel.(list{i})=1;
        elseif strfind(list{i},'fbc_version')
            sbmlModel.(list{i})= 2;
        elseif strfind(list{i},'typecode')
            sbmlModel.(list{i})='SBML_MODEL';
        elseif strfind(list{i},'fbc_activeObjective')
            sbmlModel.(list{i})='obj'; 
        elseif strfind(list{i},'id')
            sbmlModel.(list{i})='deFBAmodel';
        elseif strfind(list{i},'note')
                sbmlModel.(list{i})=['<body xmlns="http://www.w3.org/1999/xhtml"><p> </p></body>'];        
        elseif strfind(list{i},'sboTerm')
            sbmlModel.(list{i})=-1; 
        else
            sbmlModel.(list{i})=''; 
        end

    end
end

%% format ids SBML
function str = formatID(str)
str = strrep(str,'-','_DASH_');
str = strrep(str,'/','_FSLASH_');
str = strrep(str,'\','_BSLASH_');
str = strrep(str,'(','_LPAREN_');
str = strrep(str,')','_RPAREN_');
str = strrep(str,'[','_LSQBKT_');
str = strrep(str,']','_RSQBKT_');
str = strrep(str,',','_COMMA_');
str = strrep(str,'.','_PERIOD_');
str = strrep(str,'''','_APOS_');
str = regexprep(str,'\(e\)$','_e');
str = strrep(str,'&','&amp;');
str = strrep(str,'<','&lt;');
str = strrep(str,'>','&gt;');
str = strrep(str,'"','&quot;');
end