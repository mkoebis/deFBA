function [rxns] = findRxnsFromMets(model,metsToFind)
% finds IDs of reactions that produce or consume the metabolites in the
% metsToFind list
%
%INPUTS
% model         model structure containing the stoichiometric matrix as field S
%               and the metabolites field mets
% metsToFind    cell array of metabolite IDs for which the reactions should
%               be found
%
%OUTPUTS
%rxns           cell array of reaction IDs that produce or consume the metabolites in the
%               metsToFind list
%
% Alexandra Reimers 13/07/2017

    [~,idx] = ismember(metsToFind,model.mets);
    rxns = {};
    for i=1:length(idx)
        idxR = find(model.S(idx(i),:));
        rxns = [rxns;model.rxns(idxR)];
    end
    rxns = unique(rxns);
end