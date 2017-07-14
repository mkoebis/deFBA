function model = readSBML_ram(fileName,externalCompartmentID, biomassCompartmentID, storageCompartmentID)

% readSBML_adjusted reads in a SBML model in ram format as a deFBA matlab structure
%
%
%INPUTS
% fileName                  File name for file to read in
% externalCompartmentID     String containting the id of the external
%                           metabolites compartment
% biomassCompartmentID      String containting the id of the biomass compartment
% storageCompartmentID      String containting the id of the storage compartment
%
%OUTPUT
% model                     deFBA model structure with the fields:
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
%   quotaWeights                array storing the molecular weights of the quota in Da
%   rxnEnzRules                 length(rxns)xlength(enz) 0-1 matrix describing if enzyme j catalysises reaction i (entry at row i, column j is 1) or not (entry at row i, column j is 0)
%   Kcat_f                      length(rxns)xlength(enz) matrix storing the kcat of enzyme j for the forward direction of reaction i (entry at row i, column j is kcat)
%   quotaInitial                array containing the total biomass percentages that have to be satisfied for each quota compound at each time point
%   initialBiomass              initial values for the quota and enzymes amounts
%   Y0                          initial values for the storage (first noStorage entries) and external metabolites (rest) amounts
%   maintenanceID               ID of the maintenance reaction
%   maintenanceValue            maintenance dependence on total biomass
%   gprComp                     cell array that specifies for each reaction the recipe for building its corresponding enzyme from the individual gene ids ( e.g. 3*gene1 AND 2*gene2 means that the enzyme is made of three copies of gene1 and 2 copies of gene2)
%
% Alexandra Reimers 13/07/2017

%% call the TranslateSBML function from the libSBML to get a Matlab SBML structure

modelSBML = TranslateSBML(fileName); 

nMets = length(modelSBML.species);
nRxns = length(modelSBML.reaction);

%% get initial metabolite list
speciesList = {};
tmpSpecies = [];

for i = 1:nMets
    % ignore boundary metabolites that are constant
    if (~(modelSBML.species(i).boundaryCondition && modelSBML.species(i).constant))
        tmpSpecies = [tmpSpecies  modelSBML.species(i)];
        speciesList{end+1} = modelSBML.species(i).id;         
    end
end

nMets = length(speciesList);

%% build stoichiometric matrix, reaction list and rxnGene matrix
S = sparse(nMets,nRxns);
spontaneousRxn = zeros(nRxns,1);
rev = zeros(nRxns,1);
rxns = cell(nRxns,1);

%% build metabolite list and get sizes
mets = cell(nMets,1);
metNames = cell(nMets,1);
Ymet = [];
quotaMet = [];
enzMet = [];
storage = [];

for i = 1:nMets
    % parse metabolite id
    metID = speciesList{i};
    idx = strfind(metID,'_');
    if ~isempty(idx)
        idx = idx(end);
        metID = [metID(1:idx-1),'[',metID(idx+1:end),']'];
    end
    mets{i} = metID;
    metNames{i} = tmpSpecies(i).name;
     
    if strcmp(tmpSpecies(i).compartment,biomassCompartmentID)
        if ~ismember(tmpSpecies(i).id,{modelSBML.fbc_geneProduct.fbc_id})
            quotaMet = [quotaMet,i];
        else
            enzMet = [enzMet,i];
        end
    end
    if strcmp(tmpSpecies(i).compartment,externalCompartmentID)
        Ymet = [Ymet,i];
    end
    if strcmp(tmpSpecies(i).compartment,storageCompartmentID)
        storage = [storage,i];
    end
end

Xmet =  setdiff(1:1:length(tmpSpecies),[Ymet,quotaMet,enzMet,storage]);

sizeXmet = length(Xmet);
sizePmet = length(quotaMet)+length(enzMet);
sizeQuotaMet = length(quotaMet);
sizeYmet = length(storage) +length(Ymet);
noStorage = length(storage);

% get quotaInitial
quota_tmp = tmpSpecies(quotaMet);
quotaInitial = zeros(1,sizeQuotaMet);
quotaWeights = zeros(1,sizeQuotaMet);
for i=1:length(quota_tmp)
    [mw,bp] = parseSpeciesAnnotation(quota_tmp(i),modelSBML.parameter);
    quotaInitial(i) = bp;
    quotaWeights(i) = mw;
end

% get proteinWeights
prot_tmp = tmpSpecies(enzMet);
stor_tmp = tmpSpecies(storage);
proteinWeights = zeros(length(prot_tmp),1);
storageWeight = zeros(length(stor_tmp),1);
for i=1:length(prot_tmp)
    [mw,~] = parseSpeciesAnnotation(prot_tmp(i),modelSBML.parameter);
    proteinWeights(i) = mw;
end
for i=1:length(stor_tmp)
    [mw,~] = parseSpeciesAnnotation(stor_tmp(i),modelSBML.parameter);
    storageWeight(i) = mw;
end

% build enzyme list and rxn-enzyme association matrix
enz = {tmpSpecies(enzMet).id};
rxnEnzRules = zeros(nRxns,length(enz));
rxnGeneMat = [];
nGenes = 0;
ecNumbers = cell(nRxns,1);
genes = {};
grRules = cell(nRxns,1);
rxnNames = cell(nRxns,1);
gprComp = cell(nRxns,1);

for i = 1:nRxns
    %% read EC number
    notesField = modelSBML.reaction(i).notes;
    if (~isempty(notesField))
        ecNumbers{i} = extractEC(notesField);
    end

    %% read the gene product association and at the same time build rxnEnzRules,enz,genes,rxnGeneMat
    if isfield(modelSBML.reaction,'fbc_geneProductAssociation')
        if size(modelSBML.reaction(i).fbc_geneProductAssociation,2)~=0 % the Matlab structure of the "geneAssociation" is defined.
            if size(modelSBML.reaction(i).fbc_geneProductAssociation.fbc_association,2)~=0
                tmpGeneProduct = modelSBML.reaction(i).fbc_geneProductAssociation.fbc_association.fbc_association;
                gprComp{i} = tmpGeneProduct;
                if ~strcmp(tmpGeneProduct,'')
                    if isempty(strfind(tmpGeneProduct,' AND '))% name given directly by gene -> strip the stoichiometry
                        tmpName = tmpGeneProduct(strfind(tmpGeneProduct,'*')+1:end);

                        %strip potential underscores and number after that
                        if ~isempty(strfind(tmpName,'_'))
                            tmpName = tmpName(1:(strfind(tmpName,'_')-1));
                        end

                        if ~ismember(tmpName,genes)
                            genes{nGenes+1} = tmpName;
                            nGenes = nGenes+1;
                            grRules{i} = tmpName;
                            tmpUnitVector = zeros(nRxns,1);
                            tmpUnitVector(i) = 1;
                            rxnGeneMat = [rxnGeneMat,tmpUnitVector];
                        else
                            grRules{i} = tmpName;                      
                            [~,idx] = ismember(tmpName,genes);
                            rxnGeneMat(i,idx) = 1;                        
                        end
                        idx = ~cellfun(@isempty,strfind(enz,tmpName));
                        % check if it catalyses rxns in different compartments
                        if length(find(idx))>1
                            s = find(idx);
                            for h=1:length(s)
                                % assign it to the first that is not yet assigned and skip the rest
                                if sum(rxnEnzRules(:,s(h))) == 0
                                    rxnEnzRules(i,s(h)) = 1;
                                    break;
                                end
                            end
                        else
                            rxnEnzRules(i,idx) = 1;
                        end
                    else % name is a complex
                        for k=1:length(modelSBML.fbc_geneProduct)
                            if strcmp(modelSBML.fbc_geneProduct(k).fbc_label,tmpGeneProduct)
                                tmpName = modelSBML.fbc_geneProduct(k).fbc_id;
                                break;
                            end
                        end
                        idx = ~cellfun(@isempty,strfind(enz,tmpName));
                        % check if it catalyses rxns in different compartments
                        if length(find(idx))>1
                            s = find(idx);
                            for h=1:length(s)
                                % assign it to the first that is not yet assigned and skip the rest
                                if sum(rxnEnzRules(:,s(h))) == 0
                                    rxnEnzRules(i,s(h)) = 1;
                                    break;
                                end
                            end
                        else
                            rxnEnzRules(i,idx) = 1;
                        end

                        % get gene names and add them
                        idx = k;
                        gr = modelSBML.fbc_geneProduct(idx).fbc_label;
                        tokens = strsplit(gr,' AND ');
                        grRules{i} = '';

                        for z = 1:length(tokens)
                            t = tokens{z};
                            tmpNameG = t(cell2mat(strfind(tokens(z),'*'))+1:end);
                            if ~ismember(tmpNameG,genes)
                                genes{nGenes+1} = tmpNameG;
                                nGenes = nGenes+1;
                                tmpUnitVector = zeros(nRxns,1);
                                tmpUnitVector(i) = 1;
                                rxnGeneMat = [rxnGeneMat,tmpUnitVector];
                            else
                                [~,idx] = ismember(tmpNameG,genes);
                                rxnGeneMat(i,idx) = 1;
                            end
                            grRules{i} = [grRules{i},tmpNameG,' AND '];
                        end
                        tmpString = grRules{i};
                        grRules{i} = tmpString(1:end-5);
                    end
                else % if it doesn't have a gene association then it is spontaneous
                    spontaneousRxn(i) = 1;
                end
            end
        else
             spontaneousRxn(i) = 1;
        end
    end
    
    rev(i) = modelSBML.reaction(i).reversible;
    rxnNames{i} = modelSBML.reaction(i).name;
    rxns{i} = modelSBML.reaction(i).id;
    
    % Construct S-matrix
    reactantStruct = modelSBML.reaction(i).reactant;
    for j = 1:length(reactantStruct)
        speciesID = find(strcmp(reactantStruct(j).species,speciesList));
        if (~isempty(speciesID))
            stoichCoeff = reactantStruct(j).stoichiometry;
            S(speciesID,i) = -stoichCoeff;
        end
    end
    productStruct = modelSBML.reaction(i).product;
    for j = 1:length(productStruct)
        speciesID = find(strcmp(productStruct(j).species,speciesList));
        if (~isempty(speciesID))
            stoichCoeff = productStruct(j).stoichiometry;
            S(speciesID,i) = stoichCoeff;
        end
    end    
end

%% read kcat values
kcat = rxnEnzRules;
tmp_rxns = modelSBML.reaction;
for i=1:nRxns
   [maintenanceScale, kplus, kminus] = parseReactionAnnotation(tmp_rxns(i),modelSBML.parameter);
   kcat(i,logical(rxnEnzRules(i,:))) = kplus;
   
   if maintenanceScale>0 || ~isempty(strfind(tmp_rxns(i).id,'maintenance'))
       maintenanceID = tmp_rxns(i).id;
       maintenanceValue = maintenanceScale;
   end
end

%% correct enzyme names
for i=1:length(enz)
    idx = find(rxnEnzRules(:,i));
    str = 'E';
    for k=1:length(idx)
        str = [str,'_',rxns{idx(k)}];
    end
    str = [str,'[bio]'];
    
    enz{i} = str;
    mets{i+sizeXmet+sizeYmet+sizeQuotaMet} = str;
end

%% collect everything into a structure
model.rxns = rxns;
model.mets = mets;
model.S = S;
model.genes = columnVector(genes);
model.N = 10;
model.epsilon = 1e-4;
model.beta = 0;
model.sizeXmet = sizeXmet;
model.sizePmet = sizePmet;
model.sizeQuotaMet = sizeQuotaMet;
model.sizeYmet = sizeYmet;
model.noStorage = noStorage;
model.quotaInitial = quotaInitial;
model.proteinWeights = proteinWeights;
model.quotaWeights = quotaWeights;
model.storageWeight = storageWeight;
model.enz = enz;
model.rev = rev;
model.gprComp = gprComp;
model.Kcat_f = sparse(kcat);
model.rxnEnzRules = rxnEnzRules;
model.initialBiomass = [tmpSpecies([quotaMet,enzMet]).initialAmount];
model.Y0 = [tmpSpecies(storage).initialAmount  tmpSpecies(Ymet).initialAmount];
model.carbonSources = model.mets(Ymet);        
model.rxnGeneMat = rxnGeneMat;
model.spontaneousRxn = spontaneousRxn;

%% reorder S for rxns
exc=[];
for i=1:nRxns
    if (sum(full(S(:,i)>=0))==length(mets))
        exc = [exc,i];
    end
    if (sum(full(S(:,i)<=0))==length(mets))
        exc = [exc,i];
    end  
end

exc = [exc,findRxnIDs(model,findRxnsFromMets(model,model.mets([storage,Ymet])))'];
quotaRxn = findRxnIDs(model,findRxnsFromMets(model,model.mets(quotaMet)))';
enzRxn = findRxnIDs(model,findRxnsFromMets(model,model.mets(enzMet)))';
Xrxn = setdiff(1:1:nRxns,[exc,quotaRxn,enzRxn]);


%% reorder S for metabolites
model.S = model.S([Xmet,storage,setdiff(Ymet,storage),quotaMet,enzMet],:);
model.S = model.S(:,[Xrxn,exc,quotaRxn,enzRxn]);

model.sizeXrxn = length(Xrxn);
model.sizeYrxn = length(exc);
model.sizeQuotaRxn = length(quotaRxn);
model.sizePrxn = length(quotaRxn) + length(enzRxn);

model.rxns = model.rxns([Xrxn,exc,quotaRxn,enzRxn]);
model.rev = model.rev([Xrxn,exc,quotaRxn,enzRxn]);
model.gprComp = model.gprComp([Xrxn,exc,quotaRxn,enzRxn]);
model.rxnNames = columnVector(rxnNames([Xrxn,exc,quotaRxn,enzRxn]));
model.grRules = columnVector(grRules([Xrxn,exc,quotaRxn,enzRxn]));
model.rxnEnzRules = model.rxnEnzRules([Xrxn,exc,quotaRxn,enzRxn],:);
if exist('ecNumbers','var')
    model.rxnECNumbers = columnVector(ecNumbers([Xrxn,exc,quotaRxn,enzRxn]));
end
model.rxnGeneMat = model.rxnGeneMat([Xrxn,exc,quotaRxn,enzRxn],:);
model.spontaneousRxn = model.spontaneousRxn([Xrxn,exc,quotaRxn,enzRxn]);
model.Kcat_f = model.Kcat_f([Xrxn,exc,quotaRxn,enzRxn],:);

model.mets = model.mets([Xmet,storage,setdiff(Ymet,storage),quotaMet,enzMet]);
model.metNames = columnVector(metNames([Xmet,storage,setdiff(Ymet,storage),quotaMet,enzMet]));
model.noRxn = length(model.rxns);

model.tf = 1;

if exist('maintenanceID','var')
    model.maintenanceID = maintenanceID;
    model.maintenanceValue = maintenanceValue;
end

end


function vec = columnVector(vec)
    %columnVector Converts a vector to a column vector
    [n,m] = size(vec);

    if (n < m)
        vec = vec';
    end
end

function [mw, bp] = parseSpeciesAnnotation(speciesEntry,parameterEntries)
    % parse annotation to get the parameter names
    str = speciesEntry.annotation;
    idx = regexp(str,'(\"[\w\s]+\")');% extract the indices of the two fields
    mw_str = str(idx(1)+1:idx(2)-1);
    bp_str = str(idx(2)+1:end);
    
    idx = strfind(mw_str,'"');
    mw_str = mw_str(1:idx-1);
    idx = strfind(bp_str,'"');
    bp_str = bp_str(1:idx-1);
    
    % extract parameter values
    [~,idx] = ismember(mw_str,{parameterEntries.id});
    mw = parameterEntries(idx).value;
    [~,idx] = ismember(bp_str,{parameterEntries.id});
    bp = parameterEntries(idx).value;
end

function [maintenanceScale,kplus, kminus] = parseReactionAnnotation(rxnEntry,parameterEntries)
    % parse annotation to get the parameter names
    str = rxnEntry.annotation;
    idx = regexp(str,'(\"[\w\s]+\")');% extract the indices of the two fields
    
    m_str = str(idx(1)+1:idx(2)-1);    
    k1_str = str(idx(2)+1:idx(3)-1);
    k2_str = str(idx(3)+1:end);
    
    idx = strfind(m_str,'"');
    m_str = m_str(1:idx-1);
    idx = strfind(k1_str,'"');
    k1_str = k1_str(1:idx-1);
    idx = strfind(k2_str,'"');
    k2_str = k2_str(1:idx-1);
    
     % extract parameter values
    [~,idx] = ismember(m_str,{parameterEntries.id});
    maintenanceScale = parameterEntries(idx).value;
    [~,idx] = ismember(k1_str,{parameterEntries.id});
    kplus = parameterEntries(idx).value;
    [~,idx] = ismember(k2_str,{parameterEntries.id});
    kminus = parameterEntries(idx).value;
end

function [ecNumber] = extractEC(notesField)
    %parseSBMLNotesField Parse the notes field of an SBML file to extract
    %the EC number
    %
    % [ecNumber] = extractEC(notesField)

    if isempty(regexp(notesField,'html:p', 'once'))
        tag = 'p';
    else
        tag = 'html:p';
    end

    ecNumber = '';

    [~,fieldList] = regexp(notesField,['<' tag '>.*?</' tag '>'],'tokens','match');

    for i = 1:length(fieldList)
        fieldTmp = regexp(fieldList{i},['<' tag '>(.*)</' tag '>'],'tokens');
        fieldStr = fieldTmp{1}{1};
        if (regexp(fieldStr,'EC Number'))
            ecNumber = regexprep(strrep(fieldStr,'EC Number:',''),'^(\s)+','');
        end
    end
end