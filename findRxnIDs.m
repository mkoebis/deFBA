function [id] = findRxnIDs(model,rxnsToFind)
% finds IDs of reactions in the rxnsToFind list
%
%INPUTS
% model         model structure containing the stoichiometric matrix as field S
%               and the reactions field rxns
% rxnsToFind    cell array of reaction IDs for which the reactions should
%               be found
%
%OUTPUTS
%id           cell array of reaction IDs for the reactions in rxnsToFind list
%
% Alexandra Reimers 13/07/2017

    [~,id] = ismember(rxnsToFind,model.rxns);
    id = unique(id);
end