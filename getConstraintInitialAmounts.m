function [Aeq,beq] = getConstraintInitialAmounts(model)
% getConstraintInitialAmounts returns, given a deFBA model structure, 
% the equality constraint matrix and its corresponding right hand side 
% for fixing the initial conditions. If no initial amounts for
% macromolecules are provided, the constraint imposes that the
% macromolecule and storage weights at the beginning have to add up to one,
% and that the percentages of total biomass for the quota are fixed at the
% initial point.
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
% Aeq         matrix containing initial amounts equality constraint for deFBA: rows = constraints, columns = LP variables
% beq         vector containing the right hand side of the initial amounts equality constraint in Aeq, s.t. Aeq*x=beq, where x is the solution vector
% 
% Alexandra Reimers 14/07/2017

    cons = emptySolutionStruct(model);
    t=1;
	if isfield(model,'initialBiomass')
		Aeq1 = sparse(model.sizePmet,size(toVector(cons,model),2));
		for l = 1:model.sizePmet
			idx = getIndexVariable(model,'p0',l);
            Aeq1(t,idx) = 1;
            t = t+1;
		end
		beq1 = model.initialBiomass;  
        beq1(model.sizeQuotaMet+1:end) = beq1(model.sizeQuotaMet+1:end)/model.epsilon;
	else
 		Aeq1 = [];
 		beq1 = [];
    end
    t = 1;
    if isfield(model,'Y0')
        Aeq2 = sparse(model.sizeYmet,size(toVector(cons,model),2));
        for l = 1:model.sizeYmet
            idx = getIndexVariable(model,'y0',l);
            Aeq2(t,idx) = 1;
            t = t+1;
        end
        beq2 = model.Y0;     
    else
        Aeq2 = [];
        beq2 = [];
    end
    

    % the constraint below should only be active when we do not impose an
    % initial biomass
    if ~isfield(model,'initialBiomass')    
        Aeq = sparse(length(model.quotaInitial)+1,size(toVector(cons,model),2));
        beq = zeros(1,length(model.quotaInitial)+1);
        % sum(quota0)+sum(prot0)+sum(storage0)=1g
        t=1;
        
        idx = getIndexVariable(model,'p0',1:1:model.sizeQuotaMet);
        Aeq(t,idx) = model.quotaWeights;

        idx = getIndexVariable(model,'p0',model.sizeQuotaMet+1:1:model.sizePmet);
        Aeq(t,idx) = model.proteinWeights*model.epsilon;

        idx = getIndexVariable(model,'y0',1:1:model.noStorage);
        if (~isempty(idx))
            Aeq(t,idx) = model.storageWeight;
        end

        beq(1) = 1;

        % p0(quota)=quotaInitial
        t=t+1;
        for l = 1:length(model.quotaInitial)
            idx = getIndexVariable(model,'p0',l);
            % only fix the nonzero ones
            if model.quotaInitial(l)~=0
                Aeq(t,idx) = 1;
            end
            t = t+1;
        end

        beq(2:end) = model.quotaInitial';
    else
        Aeq = [];
        beq = [];
    end
    
    
    Aeq = [Aeq;Aeq1;Aeq2];
    beq = [beq,beq1,beq2];
end
