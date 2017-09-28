function model = changeGeneAssociation(model,rxnName,grRule)
% Change gene associations in a model
%
% model = changeGeneAssociation(model,rxnName,grRule,geneName,systName)
%
%INPUTS
% model             model structure
% rxnName           Name of the new reaction
% grRule            Gene-reaction rule in boolean format (and/or allowed)
%
%OUTPUT
% model             model structure with new gene reaction associations
%
% Alexandra Reimers 26/09/2017

    [isInModel,rxnID] = ismember(rxnName,model.rxns);

    if (~isInModel)
        error(['Reaction ' rxnName ' not in the model']);
    end

    if ~isfield(model,'genes')
        model.genes = {};
    end
    nGenes = length(model.genes);
    model.rules{rxnID} = '';

    model.rxnGeneMat(rxnID,:) = zeros(1,nGenes);
    % Remove extra white space
    grRule = regexprep(grRule,'\s{2,}',' ');
    grRule = regexprep(grRule,'( ','(');
    grRule = regexprep(grRule,' )',')');


    if (~isempty(grRule))
        % Remove extra white space
        grRule = regexprep(grRule,'\s{2,}',' ');
        grRule = regexprep(grRule,'( ','(');
        grRule = regexprep(grRule,' )',')');
        [genes,rule] = parseBoolean(grRule);

        for i = 1:length(genes)
            geneID = find(strcmp(model.genes,genes{i}));
            if (isempty(geneID))
                model.genes = [model.genes; genes(i)];
                nGenes = length(model.genes);

                model.rxnGeneMat(rxnID,end+1) = 1;

                rule = strrep(rule,['x(' num2str(i) ')'],['x(' num2str(nGenes) ')']);
            else
                model.rxnGeneMat(rxnID,geneID) = 1;
                rule = strrep(rule,['x(' num2str(i) ')'],['x(' num2str(geneID) ')']);
            end
        end
        model.rules{rxnID} = rule;
    end

    model.grRules{rxnID} = grRule;

    %make sure variables are column vectors
    model.rules = columnVector(model.rules);
    model.grRules = columnVector(model.grRules);
end