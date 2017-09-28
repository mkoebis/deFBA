function modelOut = removeRxns(model,rxnRemoveList)
%removeRxns Remove reactions from a model
%
%INPUTS
% model             deFBA model structure
% rxnRemoveList     Cell array of reaction names to be removed
%
%OUTPUT
% model             deFBA model w/o selected reactions
%
% Alexandra Reimers 26/09/2017


    [~,nRxns] = size(model.S);
    modelOut = model;

    % Find indices of rxns in the model
    [isValidRxn,removeInd] = ismember(rxnRemoveList,model.rxns);
    removeInd = removeInd(isValidRxn);

    % Construct vector to select rxns to be included in the models
    selectRxns = true(nRxns,1);
    selectRxns(removeInd) = false;

    % Construct new model
    if isfield(model,'description')
        modelOut.description = model.description;
    end

    modelOut.S = model.S(:,selectRxns);
    modelOut.rxns = model.rxns(selectRxns);
    modelOut.rev = model.rev(selectRxns);
    
    if (isfield(model,'genes'))
        modelOut.genes = model.genes;
        modelOut.grRules = model.grRules(selectRxns);
    end
    
    if (isfield(model,'rxnGeneMat'))
        modelOut.rxnGeneMat = model.rxnGeneMat(selectRxns,:);
    end
    
    if (isfield(model,'rules'))
        modelOut.rules = model.rules(selectRxns);
    end
    
    if (isfield(model,'rxnNames'))
        modelOut.rxnNames = model.rxnNames(selectRxns);
    end
    
    if (isfield(model, 'rxnECNumbers'))
      modelOut.rxnECNumbers = model.rxnECNumbers(selectRxns);
    end

end