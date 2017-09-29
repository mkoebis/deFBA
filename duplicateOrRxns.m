function newModel = duplicateOrRxns(model)
    [~,ors,~] = parseGeneAssociation(model);
    idx = find(ors>0);
    
    newModel = model;
    
    % for each reaction with an or combination, we add as many duplicates
    % as we have disjunctions +1
    % we delete the previous reaction which contains the or combination
    
    for i=1:length(idx)
                
        % split the gene association
        mysplit = strsplit(model.grRules{idx(i)},' OR ');
        
        % remove parentheses and add reactions
        for j=1:length(mysplit)
            mysplit{j} = strrep(mysplit{j},'(','');
            mysplit{j} = strrep(mysplit{j},')','');
            
            % make sure new gene association is not empty
            assert(~strcmp(mysplit{j},''));
            
            % add reaction to new model
            rxnID = sprintf('%s_%d',model.rxns{idx(i)},j);
            rxnName = sprintf('%s isoenzyme %d',model.rxnNames{idx(i)},j);
            newModel = addReaction(newModel,{rxnID,rxnName},...
                model.mets(model.S(:,idx(i))~=0),...
                model.S(model.S(:,idx(i))~=0,idx(i)),...
                model.rev(idx(i))==1, mysplit{j});              
        end
        i
        model.rxns(idx(i))
        % remove original reaction since now we have the replacements
        newModel = removeRxns(newModel,model.rxns(idx(i)));      
    end
end