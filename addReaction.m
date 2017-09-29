function [model,rxnIDexists] = addReaction(model,rxnName,metaboliteList,stoichCoeffList,revFlag,grRule)
%addReaction Add a reaction to the model or modify an existing reaction
%
% model = addReaction(model,rxnName,metaboliteList,stoichCoeffList,revFlag,grRule)
%
%INPUTS
% model             model structure
% rxnName           Reaction id (i.e. 'ACALD')
%                   (Note: can also be a cell array {'id','name'}
% metaboliteList    Cell array of metabolite names or alternatively the
%                   reaction formula for the reaction
% stoichCoeffList   List of stoichiometric coefficients (reactants -ve,
%                   products +ve)
%
%OPTIONAL INPUTS
% revFlag           Reversibility flag (Default = true)
% grRule            Gene-reaction rule in boolean format (and/or allowed)
%                   (Default = '');
%
%OUTPUTS
% model             COBRA model structure with new reaction
% rxnIDexists       Empty if the reaction did not exist previously, or if
%                   checkDuplicate is false. Otherwise it contains the ID
%                   of an identical reaction already present in the model.
%
% Alexandra Reimers 26/09/2017

rxnIDexists = [];

if iscell(rxnName)&&length(rxnName)>1
    rxnNameFull = rxnName{2};
    rxnName = rxnName{1};
end

% Figure out if reaction already exists
nRxns = length(model.rxns);
if (sum(strcmp(rxnName,model.rxns)) > 0)
    warning('Reaction with the same name already exists in the model');
    [~,rxnID] = ismember(rxnName,model.rxns);
    oldRxnFlag = true;
else
    rxnID = nRxns+1;
    oldRxnFlag = false;
end

% Reversibility
if (nargin < 5 || isempty(revFlag))
    if (oldRxnFlag)
        revFlag = model.rev(rxnID);
    else
        revFlag = true;
    end
end


% Missing arguments
if (nargin < 6) && (isfield(model,'grRules'))
    if (oldRxnFlag)
        grRule = model.grRules{rxnID};
    else
        grRule = '';
    end
end

nMets = length(model.mets);
Scolumn = sparse(nMets,1);

modelOrig = model;

% Update model fields
model.rxns{rxnID,1} = rxnName;
if (revFlag)
    model.rev(rxnID,1) = 1;
else
    model.rev(rxnID,1) = 0;
end

if revFlag
    model.lb(rxnID,1) = -1000;
    model.ub(rxnID,1) = 1000;
else
    model.lb(rxnID,1) = 0;
    model.ub(rxnID,1) = 1000;
end

if (isfield(model,'rxnNames'))
    if exist('rxnNameFull','var')
        model.rxnNames{rxnID,1} = rxnNameFull;
    else
        model.rxnNames{rxnID,1} = model.rxns{rxnID};
    end
end

if isfield(model,'rxnECNumbers')
    model.rxnECNumbers{rxnID,1} = '';
end

%Give warning and combine the coeffeicient if a metabolite appears more than once
[metaboliteListUnique,~,IC] = unique(metaboliteList);
if numel(metaboliteListUnique) ~= numel(metaboliteList)
    warning('Repeated mets in the formula for rxn ''%s''. Combine the stoichiometry.', rxnName)
    stoichCoeffListUnique = zeros(size(metaboliteListUnique));
    for nMetsUnique = 1:numel(metaboliteListUnique)
        stoichCoeffListUnique(nMetsUnique) = sum(stoichCoeffList(IC == nMetsUnique));
    end
    %preserve the order of metabolites:
    metOrder = [];
    for i = 1:numel(IC)
        if ~ismember(IC(i), metOrder)
            metOrder = [metOrder; IC(i)];
        end
    end
    metaboliteList = metaboliteListUnique(metOrder);
    stoichCoeffList = stoichCoeffListUnique(metOrder);
end

% Figure out which metabolites are already in the model
[isInModel,metID] = ismember(metaboliteList,model.mets);
nNewMets = sum(~isInModel);

% Construct S-matrix column

for i = 1:length(metaboliteList)
    if (isInModel(i))
        Scolumn(metID(i),1) = stoichCoeffList(i);
    else
        warning(['Metabolite ' metaboliteList{i} ' not in model - added to the model']);
        Scolumn(end+1,1) = stoichCoeffList(i);
        model.mets{end+1,1} = metaboliteList{i};

        if (isfield(model,'metNames'))      %Prompts to add missing info if desired
            model.metNames{end+1,1} = regexprep(metaboliteList{i},'(\[.+\]) | (\(.+\))','') ;
            warning(['Metabolite name for ' metaboliteList{i} ' set to ' model.metNames{end}]);
        end
    end
end

if (isfield(model,'genes'))
    if (nargin < 11)
        model = changeGeneAssociation(model,rxnName,grRule);
    else
        %fprintf('In addReaction, the class of systNameList is %s',
        %class(systNameList)); % commented out by Thierry Mondeel
        model = changeGeneAssociation(model,rxnName,grRule,geneNameList,systNameList);
    end
end

% Figure out if the new reaction already exists
rxnInModel=false;
if (nNewMets > 0) && isempty(find(newMetsCoefs == 0, 1))
    Stmp = [model.S;sparse(nNewMets,nRxns)];
else
    Stmp = model.S;
end

if (oldRxnFlag)
    model.S = Stmp;
    model.S(:,rxnID) = Scolumn;
else
    model.S = [Stmp Scolumn];
end

end
