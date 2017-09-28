function model = removeMetabolites(model,metaboliteList)
%removeMetabolites Remove metabolites from a model
%
% model = removeMetabolites(model,metaboliteList,removeRxnFlag)
%
%INPUTS
% model             deFBA model structure
% metaboliteList    List of metabolites to be removed
%
%OUTPUT
% model             deFBA model with removed metabolites
%
% Alexandra Reimers 26/09/2017

selMets = ~ismember(model.mets,metaboliteList);

model.S = model.S(selMets,:);

model.mets = model.mets(selMets);
model.metNames = model.metNames(selMets);

end