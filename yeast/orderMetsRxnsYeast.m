function model = orderMetsRxnsYeast(model)
  
    % get index of glycogen, trehalose and quota
    glycogenIdx = findMetIDs(model,'s_0773[c03]');
    trehaloseIdx = findMetIDs(model,'s_1520[c03]');
    glucoseIdx = findMetIDs(model,'s_0565[c06]');
    galactoseIdx = findMetIDs(model,'s_0559[c06]');
    ethanolIdx = findMetIDs(model,'s_0681[c06]');
    quotaIdx = findMetIDs(model,{'B1[c03]','B2[c03]','B3[c03]','B4[c03]','B5[c03]','B6[c03]'})';
    % get index of internal metabolites
    intMetIdx = (find(~ismember([1:1:length(model.mets)],[quotaIdx;glycogenIdx;trehaloseIdx;glucoseIdx;galactoseIdx;ethanolIdx])))';
    
    %% reorder metabolites
    model.mets = model.mets([intMetIdx;glycogenIdx;trehaloseIdx;glucoseIdx;galactoseIdx;ethanolIdx;quotaIdx]);
    model.metNames = model.metNames([intMetIdx;glycogenIdx;trehaloseIdx;glucoseIdx;galactoseIdx;ethanolIdx;quotaIdx]);
    model.S = model.S([intMetIdx;glycogenIdx;trehaloseIdx;glucoseIdx;galactoseIdx;ethanolIdx;quotaIdx],:);
    model.noStorage = 2;
    %% get reaction indices
    % find transport reactions indices and quota forming reactions
    trRxn = findRxnIDs(model,model.rxns(findExcRxns(model)));
    % plus the glycogen, trehalose and ethanol forming and consuming reactions
    trRxn = [trRxn;(findRxnIDs(model,findRxnsFromMets(model,'s_0773[c03]')));...
        (findRxnIDs(model,findRxnsFromMets(model,'s_1520[c03]')));
        (findRxnIDs(model,findRxnsFromMets(model,'s_0565[c06]')));
        (findRxnIDs(model,findRxnsFromMets(model,'s_0559[c06]')));
        (findRxnIDs(model,findRxnsFromMets(model,'s_0681[c06]')))];
    
    poolRxn = [];
    c = quotaIdx;
    for i=1:length(c)
        poolRxn = [poolRxn,find(model.S(c(i),:)>0)];
    end
    trRxn = unique(trRxn);
    poolRxn = unique(poolRxn);
    intRxnIdx = find(~ismember([1:1:length(model.rxnNames)],[poolRxn';trRxn]));
    
    
    % check for duplicate reactions
    for i=1:length(model.rxns)
        assert(length(find(ismember(model.rxns,model.rxns(i))))==1)
    end
    for i=1:length(trRxn)
        assert(length(find(ismember(trRxn,trRxn(i))))==1)
    end
    for i=1:length(poolRxn)
        assert(length(find(ismember(poolRxn,poolRxn(i))))==1)
    end
    for i=1:length(intRxnIdx)
        assert(length(find(ismember(intRxnIdx,intRxnIdx(i))))==1)
    end
    
    %% reorder reactions
    model.rxns = [model.rxns(intRxnIdx);model.rxns(trRxn);model.rxns(poolRxn)];
    model.rxnNames = [model.rxnNames(intRxnIdx);model.rxnNames(trRxn);model.rxnNames(poolRxn)];
    model.S = [model.S(:,intRxnIdx),model.S(:,trRxn),model.S(:,poolRxn)];
    model.rxnGeneMat = [model.rxnGeneMat(intRxnIdx,:);model.rxnGeneMat(trRxn,:);model.rxnGeneMat(poolRxn,:)];
    model.rev = [model.rev(intRxnIdx);model.rev(trRxn);model.rev(poolRxn)];
    model.rules = [model.rules(intRxnIdx);model.rules(trRxn);model.rules(poolRxn)];
    model.grRules = [model.grRules(intRxnIdx);model.grRules(trRxn);model.grRules(poolRxn)];
%     model.rxnECNumbers = [model.rxnECNumbers(intRxnIdx);model.rxnECNumbers(trRxn);model.rxnECNumbers(poolRxn)];
    model.genes = model.genes;
    
    % adjust sizes
    model.sizeXmet = length(intMetIdx);
    model.sizeXrxn = length(intRxnIdx);
    model.sizeYmet = 5;
    model.sizeYrxn = length(trRxn);
    model.sizeQuotaMet = length(c);
    model.sizeQuotaRxn = length(poolRxn);
end
