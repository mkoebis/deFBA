function model = readSBML_ram(fileName)

% readSBML_adjusted reads in a SBML model in ram format as a deFBA matlab structure
%
%
%INPUTS
% fileName                  File name for file to read in
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
xmlModel = getAnnotations(fileName);

nMets = length(modelSBML.species);
nRxns = length(modelSBML.reaction);

%% get initial metabolite list
speciesList = {};
metCompartments = {};
tmpSpecies = [];

for i = 1:nMets
    % ignore boundary metabolites that are constant
    if (~(modelSBML.species(i).boundaryCondition && modelSBML.species(i).constant))
        tmpSpecies = [tmpSpecies  modelSBML.species(i)];
        speciesList{end+1} = modelSBML.species(i).id;  
        metCompartments{end+1} = modelSBML.species(i).compartment;
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
Xmet = [];
quotaMet = [];
enzMet = [];
storage = [];

for i = 1:nMets
    mets{i} = speciesList{i};
    metNames{i} = tmpSpecies(i).name;
      
    if ~isempty(xmlModel.mets(i).annotation)
        for j=1:length(xmlModel.mets(i).annotation)
            if strcmp(xmlModel.mets(i).annotation(j).Name,'ram:speciesType')
                if strcmp(xmlModel.mets(i).annotation(j).Value,'enzyme')
                    enzMet = [enzMet,i];
                elseif strcmp(xmlModel.mets(i).annotation(j).Value,'storage')
                    storage = [storage,i];
                elseif strcmp(xmlModel.mets(i).annotation(j).Value,'quota')
                    % check for mislabeling
                    if ismember(xmlModel.mets(i).id,{xmlModel.geneProd.associatedSpecies})
                        enzMet = [enzMet,i];
                    else                    
                        quotaMet = [quotaMet,i];
                    end
                elseif strcmp(xmlModel.mets(i).annotation(j).Value,'metabolite')
                    Xmet = [Xmet,i];
                elseif strcmp(xmlModel.mets(i).annotation(j).Value,'extracellular')
                    if (tmpSpecies(i).boundaryCondition==0)
                        Ymet = [Ymet,i];
                    end
                end
            end
        end
    end
end

sizeXmet = length(Xmet);
sizePmet = length(quotaMet)+length(enzMet);
sizeQuotaMet = length(quotaMet);
sizeYmet = length(storage) +length(Ymet);
noStorage = length(storage);

% get quotaInitial
quota_tmp = [xmlModel.mets([quotaMet,enzMet])];
quotaInitial = zeros(1,sizePmet);
quotaWeights = zeros(1,sizeQuotaMet);
objectiveWeights = zeros(1,sizePmet);
c = 0;
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
                    if ismember(quota_tmp(i).id,mets(quotaMet))
                        c = c+1;
                        quotaWeights(c) = str2double(xmlModel.params(k).value);
                    end
                end
            end
        end
        if strcmp(quota_tmp(i).annotation(j).Name,'ram:objectiveWeight')
            for k=1:length(xmlModel.params)
                if strcmp(xmlModel.params(k).id,quota_tmp(i).annotation(j).Value)
                    objectiveWeights(i) = str2double(xmlModel.params(k).value);
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
fbcGeneIDs = {xmlModel.geneProd.id};
for i = 1:nRxns
    %% read EC number
    if isfield(modelSBML.reaction(i),'cvterms')
        if ~isempty(modelSBML.reaction(i).cvterms)
            ecNumbers{i} = parseIdentifiersLink(modelSBML.reaction(i).cvterms.resources);
        end
    end
    
    %% read the gene product association and at the same time build rxnEnzRules,enz,genes,rxnGeneMat    
    if ~isempty(xmlModel.rxns(i).geneProductAssociation) 
        % get associated species
        idx = ismember(fbcGeneIDs,xmlModel.rxns(i).geneProductAssociation);
        e = xmlModel.geneProd(idx).associatedSpecies;
        gprComp{i} = xmlModel.geneProd(idx).label;
        if isempty(strfind(xmlModel.geneProd(idx).label,' AND '))% only one gene involved
            tmpName = xmlModel.geneProd(idx).label(strfind(xmlModel.geneProd(idx).label,'*')+1:end);

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
            
            [~,idx] = ismember(enz,e);            
            rxnEnzRules(i,logical(idx)) = 1;
        else % complex
            %split label
            g = strsplit(xmlModel.geneProd(idx).label,' AND ');
            
            for k=1:length(g)
                tmpName = g{k};
                tmpName = tmpName(strfind(tmpName,'*')+1:end);
                g{k} = tmpName;
                
                if ~ismember(tmpName,genes)
                    genes{nGenes+1} = tmpName;
                    nGenes = nGenes+1;
                    %grRules{i} = tmpName;
                    tmpUnitVector = zeros(nRxns,1);
                    tmpUnitVector(i) = 1;
                    rxnGeneMat = [rxnGeneMat,tmpUnitVector];
                else
                    %grRules{i} = tmpName;                      
                    [~,idx] = ismember(tmpName,genes);
                    rxnGeneMat(i,idx) = 1;                        
                end
            end
            
            tmpGrRule = strjoin(g, ' AND ');
            grRules{i} = tmpGrRule;
            
            [~,idx] = ismember(enz,e);            
            rxnEnzRules(i,logical(idx)) = 1;           
        end
    else
         spontaneousRxn(i) = 1;
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

%% read kcat values and maintenance
kcat = rxnEnzRules;
tmp_rxns = modelSBML.reaction;
for i=1:nRxns
    ann = {xmlModel.rxns(i).annotation.Name};
    idx = ismember(ann,'ram:kcatForward');
    p = ismember({xmlModel.params.id},xmlModel.rxns(i).annotation(idx).Value);
    if ~isempty(find(p))
        kplus = str2double(xmlModel.params(p).value);
    else
        kplus = 0;
    end
    
    
    idx = ismember(ann,'ram:kcatBackward');
    p = ismember({xmlModel.params.id},xmlModel.rxns(i).annotation(idx).Value);
    if ~isempty(find(p))
        kminus = str2double(xmlModel.params(p).value);
    else
        kminus = 0;
    end
    
    idx = ismember(ann,'ram:maintenanceScaling');
    p = ismember({xmlModel.params.id},xmlModel.rxns(i).annotation(idx).Value);
    maintenanceScale = str2double(xmlModel.params(p).value);
       
   kcat(i,logical(rxnEnzRules(i,:))) = kplus;
   
   if maintenanceScale>0 || ~isempty(strfind(tmp_rxns(i).id,'maintenance'))
       maintenanceID = tmp_rxns(i).id;
       maintenanceValue = maintenanceScale;
   end
end

%% collect everything into a structure
model.rxns = rxns;
model.mets = mets;
model.metCompartments = metCompartments;
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
model.objectiveWeights = objectiveWeights;
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
model.metCompartments = model.metCompartments([Xmet,storage,setdiff(Ymet,storage),quotaMet,enzMet]);
model.metNames = columnVector(metNames([Xmet,storage,setdiff(Ymet,storage),quotaMet,enzMet]));
model.noRxn = length(model.rxns);


model.tf = 1;

if exist('maintenanceID','var')
    model.maintenanceID = maintenanceID;
    model.maintenanceValue = maintenanceValue;
end

warning('Please set the ribosome ID in the field model.ribosomeID!')
warning('If there is *noncatalytic* protein quota, please add the ID of the reaction that produces it in model.ProtQuotaProdID, else add an empty string in that field')


end


function vec = columnVector(vec)
    %columnVector Converts a vector to a column vector
    [n,m] = size(vec);

    if (n < m)
        vec = vec';
    end
end

function ec = parseIdentifiersLink(link)
    ec = regexp(link{1},'http://identifiers\.org/ec-code/(?<ec>.+)','names');
    ec = ec.ec;
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
                
        
        chld = {xmlModel.rxns(i).Children.Name};
        [~,idx] = ismember('fbc:geneProductAssociation',chld);
        if idx~=0
            xmlModel.rxns(i).geneProductAssociation = xmlModel.rxns(i).Children(idx);
        else
            xmlModel.rxns(i).geneProductAssociation = '';
        end
        
        if ~isempty(xmlModel.rxns(i).geneProductAssociation)
            chld = {xmlModel.rxns(i).geneProductAssociation.Children.Name};
            [~,idx] = ismember('fbc:geneProductRef',chld);
            xmlModel.rxns(i).geneProductAssociation = xmlModel.rxns(i).geneProductAssociation.Children(idx).Attributes.Value;
        end
               
        xmlModel.rxns(i).Children = xmlModel.rxns(i).Children(2:2:end);
    end
    
    for i=1:length(xmlModel.geneProd)
        att = {xmlModel.geneProd(i).Attributes.Name};
        idx = ismember(att, 'fbc:associatedSpecies');
        xmlModel.geneProd(i).associatedSpecies = xmlModel.geneProd(i).Attributes(idx).Value;
        idx = ismember(att, 'fbc:id');
        xmlModel.geneProd(i).id = xmlModel.geneProd(i).Attributes(idx).Value;
        idx = ismember(att, 'fbc:label');
        xmlModel.geneProd(i).label = xmlModel.geneProd(i).Attributes(idx).Value;
    end
    
    xmlModel.rxns = rmfield(xmlModel.rxns,'Name');
    xmlModel.rxns = rmfield(xmlModel.rxns,'Data');
    xmlModel.rxns = rmfield(xmlModel.rxns,'Children');
end
