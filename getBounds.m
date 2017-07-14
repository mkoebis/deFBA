function [lb,ub] = getBounds(model)
% getBounds returns the lower and upper variable bounds needed to run deFBA given a deFBA model structure
%
%INPUTS
% model             deFBA model structure with the fields:
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
%   rxnEnzRules                 length(rxns)xlength(enz) 0-1 matrix describing if enzyme j catalysises reaction i (entry at row i, column j is 1) or not (entry at row i, column j is 0)
%   Kcat_f                      length(rxns)xlength(enz) matrix storing the kcat of enzyme j for the forward direction of reaction i (entry at row i, column j is kcat)
%   quotaInitial                array containing the total biomass percentages that have to be satisfied for each quota compound at each time point
%   initialBiomass              initial values for the quota and enzymes amounts
%   Y0                          initial values for the storage (first noStorage entries) and external metabolites (rest) amounts
%   maintenanceID               ID of the maintenance reaction
%   maintenanceValue            maintenance dependence on total biomass
%   gprComp                     cell array that specifies for each reaction the recipe for building its corresponding enzyme from the individual gene ids ( e.g. 3*gene1 AND 2*gene2 means that the enzyme is made of three copies of gene1 and 2 copies of gene2)
%
%OUTPUT
% lb            vector containing all variable lower bounds
% ub            vector containing all variable upper bounds
% 
% Alexandra Reimers 14/07/2017


    cons = emptySolutionStruct(model);
    lb = -inf*ones(1,size(toVector(cons,model),2));
    ub = inf*ones(1,size(toVector(cons,model),2));
    assert(length(toVector(emptySolutionStruct(model),model)) == length(lb))
    assert(length(lb) == length(ub))
    
    % positive concentrations y,p,p0,y0
    lb(getIndexVariable(model,'y0',1:1:model.sizeYmet)) = 0;
    lb(getIndexVariable(model,'p0',1:1:model.sizePmet)) = 0;
    lb(getIndexVariable(model,'y',1:1:model.N,1:1:model.sizeYmet)) = 0;
    lb(getIndexVariable(model,'p',1:1:model.N,1:1:model.sizePmet)) = 0;
    
    % helper variables for reversible reactions should be positive
    lb(getIndexVariable(model,'vRev',1:1:model.N,1:1:sum(model.rev))) = 0;
    
    % reaction irreversibility
    for i=1:model.N
        for k=1:model.noRxn
            if (model.rev(k)==0)
                lb(getIndexVariable(model,'v',i,k)) = 0;
            end
        end
    end
          
    % impose additional reaction bounds if they exist
    if isfield(model,'rxnExtraBounds')
        % impose customize lower bounds for uptake  
        for i=1:length(model.rxnExtraBounds)
            idxTime = find(model.rxnExtraBounds{i}.timePoint<=model.N);
            if ~isempty(idxTime)
                idx = getIndexVariable(model,'v',model.rxnExtraBounds{i}.timePoint(idxTime),findRxnIDs(model,model.rxnExtraBounds{i}.rxnID));
                lb(idx) = model.rxnExtraBounds{i}.lb;
                ub(idx) = model.rxnExtraBounds{i}.ub;
            end
        end
    end
        
    % impose switches in nutrients if they exist
    if isfield(model,'switches')
        for i=1:length(model.switches)
            if model.switches(i).n<=model.N
                lb(getIndexVariable(model,'y',model.switches(i).n,model.switches(i).idx)) = model.switches(i).amount;
                ub(getIndexVariable(model,'y',model.switches(i).n,model.switches(i).idx)) = model.switches(i).amount;
            end
        end
    end
end