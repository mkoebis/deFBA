function model = readSBML_ram(fileName,externalCompartmentID,renameEnzymesFlag)

% readSBML_adjusted reads in a SBML model in ram format as a deFBA matlab structure
%
%
%INPUTS
% fileName                  File name for file to read in
% externalCompartmentID     String containting the id of the external
%                           metabolites compartment
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
xmlModel = getAnnotations(fileName);

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
cmps = cell(nMets,1);
metNames = cell(nMets,1);
Ymet = [];
quotaMet = [];
enzMet = [];
storage = [];

for i = 1:nMets
    % parse metabolite id
    metID = speciesList{i};

    [m,c] = parseMetID(formatID(metID));
    
    if ~strcmp(m,'[]')&&~strcmp(c,'[]')
        mets{i} = m;
        cmps{i} = c;
    else
        mets{i} = metID;
        cmps{i} = '';
    end
    metNames{i} = tmpSpecies(i).name;
    
    if strcmp(tmpSpecies(i).compartment,externalCompartmentID)
        if (tmpSpecies(i).boundaryCondition==1)
            Ymet = [Ymet,i];
        end
    end
    
    if ~isempty(xmlModel.mets(i).annotation)
        for j=1:length(xmlModel.mets(i).annotation)
            if strcmp(xmlModel.mets(i).annotation(j).Name,'ram:speciesType')
                if strcmp(xmlModel.mets(i).annotation(j).Value,'enzyme')
                    enzMet = [enzMet,i];
                elseif strcmp(xmlModel.mets(i).annotation(j).Value,'storage')
                    storage = [storage,i];
                elseif strcmp(xmlModel.mets(i).annotation(j).Value,'quota')
                    quotaMet = [quotaMet,i];
                end
            end
        end
    end
end

Xmet =  setdiff(1:1:length(tmpSpecies),[Ymet,quotaMet,enzMet,storage]);

sizeXmet = length(Xmet);
sizePmet = length(quotaMet)+length(enzMet);
sizeQuotaMet = length(quotaMet);
sizeYmet = length(storage) +length(Ymet);
noStorage = length(storage);

% get quotaInitial
quota_tmp = xmlModel.mets(quotaMet);
quotaInitial = zeros(1,sizeQuotaMet);
quotaWeights = zeros(1,sizeQuotaMet);
for i=1:length(quota_tmp)
    for j=1:length(quota_tmp(i).annotation)
        if strcmp(quota_tmp(i).annotation(j).Name,'ram:biomassPercentage')
            for k=1:length(xmlModel.params)
                if strcmp(xmlModel.params(k).id,quota_tmp(i).annotation(j).Value)
                    quotaInitial(i) = str2double(xmlModel.params(k).value);
                end
            end
        end
        if strcmp(quota_tmp(i).annotation(j).Name,'ram:molecularWeight')
            for k=1:length(xmlModel.params)
                if strcmp(xmlModel.params(k).id,quota_tmp(i).annotation(j).Value)
                    quotaWeights(i) = str2double(xmlModel.params(k).value);
                end
            end
        end
    end
end

% get proteinWeights
prot_tmp = xmlModel.mets(enzMet);
stor_tmp = xmlModel.mets(storage);
proteinWeights = zeros(length(prot_tmp),1);
storageWeight = zeros(length(stor_tmp),1);
for i=1:length(prot_tmp)
    for j=1:length(prot_tmp(i).annotation)
        if strcmp(prot_tmp(i).annotation(j).Name,'ram:molecularWeight')
            for k=1:length(xmlModel.params)
                if strcmp(xmlModel.params(k).id,prot_tmp(i).annotation(j).Value)
                    proteinWeights(i) = str2double(xmlModel.params(k).value);
                end
            end
        end
    end
end
for i=1:length(stor_tmp)
    for j=1:length(stor_tmp(i).annotation)
        if strcmp(stor_tmp(i).annotation(j).Name,'ram:molecularWeight')
            for k=1:length(xmlModel.params)
                if strcmp(xmlModel.params(k).id,stor_tmp(i).annotation(j).Value)
                    storageWeight(i) = str2double(xmlModel.params(k).value);
                end
            end
        end
    end
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
    if isfield(xmlModel.rxns(i),'ec')
        ecNumbers{i} = xmlModel.rxns(i).ec;
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
    rxns{i} = formatID(modelSBML.reaction(i).id);
    
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

%% read kcat values and maintenance
kcat = rxnEnzRules;
tmp_rxns = modelSBML.reaction;
for i=1:nRxns
    for j=1:length(xmlModel.rxns(i).annotation)
        if strcmp(xmlModel.rxns(i).annotation(j).Name,'ram:kcatForward')
            for k=1:length(xmlModel.params)
                if strcmp(xmlModel.params(k).id,xmlModel.rxns(i).annotation(j).Value)
                    kplus = str2double(xmlModel.params(k).value);
                end
            end
        end
        if strcmp(xmlModel.rxns(i).annotation(j).Name,'ram:kcatBackward')
            for k=1:length(xmlModel.params)
                if strcmp(xmlModel.params(k).id,xmlModel.rxns(i).annotation(j).Value)
                    kminus = str2double(xmlModel.params(k).value);
                end
            end
        end
        if strcmp(xmlModel.rxns(i).annotation(j).Name,'ram:maintenanceScaling')
            for k=1:length(xmlModel.params)
                if strcmp(xmlModel.params(k).id,xmlModel.rxns(i).annotation(j).Value)
                    maintenanceScale = str2double(xmlModel.params(k).value);
                end
            end
        end
    end
    
   kcat(i,logical(rxnEnzRules(i,:))) = kplus;
   
   if maintenanceScale>0 || ~isempty(strfind(tmp_rxns(i).id,'maintenance'))
       maintenanceID = tmp_rxns(i).id;
       maintenanceValue = maintenanceScale;
   end
end

%% correct enzyme names
if renameEnzymesFlag
    for i=1:length(enz)
        idx = find(rxnEnzRules(:,i));
        str = 'E';
        for k=1:length(idx)
            str = [str,'_',rxns{idx(k)}];
        end
        str = [str,cmps{i+sizeXmet+sizeYmet+sizeQuotaMet}];

        enz{i} = str;
        mets{i+sizeXmet+sizeYmet+sizeQuotaMet} = str;
    end
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

function xmlModel = getAnnotations(filename)
    xmlModel = parseXML(filename);
    ch = {xmlModel.Children.Name};       
    xmlModel = xmlModel.Children(ismember(ch,'model'));
    
    ch = {xmlModel.Children.Name};
    xmlModel.mets = xmlModel.Children(ismember(ch,'listOfSpecies'));
    xmlModel.rxns = xmlModel.Children(ismember(ch,'listOfReactions'));
    xmlModel.params = xmlModel.Children(ismember(ch,'listOfParameters'));
    xmlModel.geneProd = xmlModel.Children(ismember(ch,'fbc:listOfGeneProducts'));
    
    ch = {xmlModel.mets.Children.Name};
    xmlModel.mets = xmlModel.mets.Children(ismember(ch,'species'));
    
    ch = {xmlModel.rxns.Children.Name};
    xmlModel.rxns = xmlModel.rxns.Children(ismember(ch,'reaction'));
    
    ch = {xmlModel.params.Children.Name};
    xmlModel.params = xmlModel.params.Children(ismember(ch,'parameter'));
    
    ch = {xmlModel.geneProd.Children.Name};
    xmlModel.geneProd = xmlModel.geneProd.Children(ismember(ch,'fbc:geneProduct'));
    
    extmets = [];
    for i=1:length(xmlModel.mets)
        atts = {xmlModel.mets(i).Attributes.Name};
        [~,idx] = ismember('id',atts);
        xmlModel.mets(i).id = xmlModel.mets(i).Attributes(idx).Value;
        
        [~,idx] = ismember('boundaryCondition',atts);
        if strcmp(xmlModel.mets(i).Attributes(idx).Value,'true')
            [~,idx] = ismember('constant',atts);
            if strcmp(xmlModel.mets(i).Attributes(idx).Value,'true')
                extmets = [extmets,i];
            end
        end
        
        if ~isempty(xmlModel.mets(i).Children)
            xmlModel.mets(i).annotation = xmlModel.mets(i).Children(2).Children(2).Children(2).Attributes;
        end
    end
    xmlModel.mets = xmlModel.mets(setdiff(1:1:length(xmlModel.mets),extmets));
    xmlModel.mets = rmfield(xmlModel.mets,'Name');
    xmlModel.mets = rmfield(xmlModel.mets,'Data');
    xmlModel.mets = rmfield(xmlModel.mets,'Children');
    
    for i=1:length(xmlModel.params)        
        atts = {xmlModel.params(i).Attributes.Name};
        [~,idx] = ismember('id',atts);
        xmlModel.params(i).id = xmlModel.params(i).Attributes(idx).Value;

        [~,idx] = ismember('value',atts);
        xmlModel.params(i).value = xmlModel.params(i).Attributes(idx).Value;
    end
    xmlModel.params = rmfield(xmlModel.params,'Name');
    xmlModel.params = rmfield(xmlModel.params,'Data');
    xmlModel.params = rmfield(xmlModel.params,'Children');
    
    for i=1:length(xmlModel.rxns)
        atts = {xmlModel.rxns(i).Attributes.Name};
        [~,idx] = ismember('id',atts);
        xmlModel.rxns(i).id = xmlModel.rxns(i).Attributes(idx).Value;
        
        chld = {xmlModel.rxns(i).Children.Name};
        [~,idx] = ismember('annotation',chld);
        xmlModel.rxns(i).annotation = xmlModel.rxns(i).Children(idx);
        
        chld = {xmlModel.rxns(i).annotation.Children.Name};
        xmlModel.rxns(i).annotation = xmlModel.rxns(i).annotation.Children(ismember(chld,'ram:RAM'));
        
        chld = {xmlModel.rxns(i).annotation.Children.Name};
        xmlModel.rxns(i).annotation = xmlModel.rxns(i).annotation.Children(ismember(chld,'ram:reaction')).Attributes;
        
        [idx] = ismember(chld,'notes');
        if any(idx~=0)
            ec = xmlModel.rxns(i).Children(idx).Children(2).Children(2).Children.Data;
            xmlModel.rxns(i).ec = regexprep(strrep(ec,'EC Number:',''),'^(\s)+','');
        end
        
        xmlModel.rxns(i).Children = xmlModel.rxns(i).Children(2:2:end);
    end
    xmlModel.rxns = rmfield(xmlModel.rxns,'Name');
    xmlModel.rxns = rmfield(xmlModel.rxns,'Data');
    xmlModel.rxns = rmfield(xmlModel.rxns,'Children');
end

function [idsMet,comps] = parseMetID(id)
    idsMet = '';
    comps = '';
    [~,aux] = regexp(id,'(?<met>.+)_(?<comp3>.+)_(?<comp2>.+)_(?<comp1>.+)','tokens','names');
     
    if isempty(aux)
       [~,aux] = regexp(id,'(?<met>.+)_(?<comp2>.+)_(?<comp1>.+)','tokens','names');
    else
        if isnan(str2double(aux.comp3))&&isnan(str2double(aux.comp2))
            idsMet = [aux.met,'[',aux.comp3,']_[',aux.comp2,']_[',aux.comp1,']'];
            comps = ['[',aux.comp3,']_[',aux.comp2,']_[',aux.comp1,']'];
        elseif ~isnan(str2double(aux.comp3))&&isnan(str2double(aux.comp2))
            idsMet = [aux.met,'_',aux.comp3,'[',aux.comp2,']_[',aux.comp1,']'];
            comps = ['[',aux.comp2,']_[',aux.comp1,']'];
        end
    end
    
    if isempty(aux)
        [~,aux] = regexp(id,'(?<met>.+)_(?<comp1>.+)','tokens','names');
    else
        if isempty(idsMet)
            if isnan(str2double(aux.comp2))
                idsMet = [aux.met,'[',aux.comp2,']_[',aux.comp1,']'];
                comps = ['[',aux.comp2,']_[',aux.comp1,']'];
            else
                idsMet = [aux.met,'_',aux.comp2,'[',aux.comp1,']'];
                comps = ['[',aux.comp1,']'];
            end
        end
    end   
    
    if isempty(idsMet)
        idsMet = [aux.met,'[',aux.comp1,']'];
        comps = ['[',aux.comp1,']'];
    end    
end

function str = formatID(str)
str = strrep(str,'_DASH_','-');
str = strrep(str,'_FSLASH_','/');
str = strrep(str,'_BSLASH_','\');
str = strrep(str,'_LPAREN_','(');
str = strrep(str,'_RPAREN_',')');
str = strrep(str,'_LSQBKT_','[');
str = strrep(str,'_RSQBKT_',']');
str = strrep(str,'_COMMA_',',');
str = strrep(str,'_PERIOD_','.');
str = strrep(str,'_APOS_','''');
str = regexprep(str,'_e','\(e\)$');
str = strrep(str,'&amp;','&');
str = strrep(str,'&lt;','<');
str = strrep(str,'&gt;','>');
str = strrep(str,'&quot;','"');
end
