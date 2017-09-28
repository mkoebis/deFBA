function [excRxn] = findExcRxns(model)
%findExcRxns Find exchange and uptake rxns
%
% [selExc,selUpt] = findExcRxns(model,inclObjFlag,irrevFlag)
%
%INPUT
% model             model structure
%
%
%OUTPUTS
% excRxn            cell array of exchange reactions
%
% Exchange reactions only have only nonnegative or only nonpositive entries
% in their corresponding S column
%
% Alexandra Reimers 26/09/2017

idx = false(length(model.rxns),1);
for i=1:length(model.rxns)
    if (sum(full(model.S(:,i)>=0))==length(model.mets))
        idx(i) = true;
    end
    if (sum(full(model.S(:,i)<=0))==length(model.mets))
        idx(i) = true;
    end    
end
excRxn = idx;

end
    