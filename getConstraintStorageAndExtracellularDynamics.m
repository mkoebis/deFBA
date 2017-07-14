function [Aeq,beq] = getConstraintStorageAndExtracellularDynamics(model)
% getConstraintStorageAndExtracellularDynamics returns, given a deFBA model
% structure, the equality constraint matrix and its corresponding right 
% hand side for imposing that the derivatives of the storage metabolites
% and extracellular metabolites amounts are equal to the product of S and 
% their corresponding production/consumption fluxes.
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
% Aeq         matrix containing storage and extracellular metabolite dynamics equality constraint for deFBA: rows = constraints, columns = LP variables
% beq         vector containing the right hand side of the storage and extracellular metabolite dynamics equality constraint in Aeq, s.t. Aeq*x=beq, where x is the solution vector
% 
% Alexandra Reimers 14/07/2017

    sizeXmet = size(model.S,1)-model.sizeYmet-model.sizePmet;
    sizeXrxn = size(model.S,2)-model.sizeYrxn-model.sizePrxn;
    assert(sizeXmet == model.sizeXmet)
    assert(sizeXrxn == model.sizeXrxn) 
    
    cons = emptySolutionStruct(model);
    Aeq = sparse(model.N*model.sizeYmet,size(toVector(cons,model),2));
    beq = zeros(1,model.N*model.sizeYmet);
    t = 1;
    for i = 1:model.N
        for l = 1:model.sizeYmet
            % compute constraint coefficient for ydot
            idx = getIndexVariable(model,'ydot',i,l);
            Aeq(t,idx) = 1;
            % compute constraint coefficient for v
            idx = getIndexVariable(model,'v',i,(sizeXrxn+1):(sizeXrxn+model.sizeYrxn));
            Aeq(t,idx) = -model.S(sizeXmet+l,(sizeXrxn+1):(sizeXrxn+model.sizeYrxn));

            if i==1
                idx = getIndexVariable(model,'y0',l);
                Aeq(t,idx) = model.beta/2;
                idx = getIndexVariable(model,'y',i,l);
                Aeq(t,idx) = model.beta/2;
            elseif i>1
                idx = getIndexVariable(model,'y',i,l);
                Aeq(t,idx) = model.beta/2;
                idx = getIndexVariable(model,'y',i-1,l);
                Aeq(t,idx) = model.beta/2;
            end

            t = t+1;
        end
    end    
end