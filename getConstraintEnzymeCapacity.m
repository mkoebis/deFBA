function [Aineq,bineq] = getConstraintEnzymeCapacity(model)
% getConstraintEnzymeCapacity returns, given a deFBA model structure, 
% the inequality constraint matrix and its corresponding right hand side 
% for imposing that at each time point each reaction flux is bound by its
% catalyzing enzyme and the kcat
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
% Aineq         matrix containing enzyme capacity inequality constraint for deFBA: rows = constraints, columns = LP variables
% bineq         vector containing the right hand side of the enzyme capacity inequality constraint in Aineq, s.t. Aineq*x<=bineq, where x is the solution vector
% 
% Alexandra Reimers 14/07/2017

    model.h = model.tf/model.N;
    sizeXYQ = size(model.S,2)-model.sizePrxn + model.sizeQuotaRxn;
    assert(sizeXYQ == model.sizeXrxn + model.sizeYrxn + model.sizeQuotaRxn)
    revRxn = model.rxns(logical(model.rev));
    
    cons = emptySolutionStruct(model);
    Aineq = sparse(model.N*length(model.enz),size(toVector(cons,model),2));
    bineq = zeros(1,model.N*length(model.enz));
    Aineq2 = sparse(model.N*2*length(revRxn),size(toVector(cons,model),2));
    bineq2 = zeros(1,model.N*2*length(revRxn));
    d = 1;
    t = 1;
    for i = 1:model.N
        for  l=1:length(model.enz)
            % find all rxns calalysed by this enzyme
            idxRxns = find(model.rxnEnzRules(:,l));

            % get enzyme indices
            if i==1
                idxp = getIndexVariable(model,'p0',l+model.sizeQuotaMet);
            else
                idxp = getIndexVariable(model,'p',i-1,l+model.sizeQuotaMet);
            end
            idxpdot = getIndexVariable(model,'pdot',i,l+model.sizeQuotaMet);
            Aineq(t,idxp) = -1;
            Aineq(t,idxpdot) = -model.h/2;


            % special case for ribosome because of scaling
            if strcmp(model.enz(l),model.ribosomeID)
                for j=1:length(idxRxns)
                    idxv = getIndexVariable(model,'v',i,idxRxns(j));
                    % quota is not scaled by epsilon
                    if strcmp(model.rxns(idxRxns(j)),model.ProtQuotaProdID)
                        Aineq(t,idxv) = 1/(model.epsilon*model.Kcat_f(idxRxns(j),end));
                    else
                        Aineq(t,idxv) = 1/model.Kcat_f(idxRxns(j),end);
                    end
                end
            else
                for j=1:length(idxRxns)
                    %if we have a reversible reaction then we have to
                    %use the helpers
                    % assumption: kcatforward = kcatreverse
                    if model.rev(idxRxns(j))==1
                        idxHelp = find(ismember(revRxn,model.rxns(idxRxns(j))));
                        idxHelp = getIndexVariable(model, 'vRev', i,idxHelp);
                        Aineq(t,idxHelp) = 1/(model.Kcat_f(idxRxns(j),l)*model.epsilon);

                        % make sure helpers bound absolute value of v
                        Aineq2(d,idxHelp) = -1;
                        Aineq2(d,getIndexVariable(model,'v',i,idxRxns(j))) = -1;
                        d = d+1;
                        Aineq2(d,idxHelp) = -1;
                        Aineq2(d,getIndexVariable(model,'v',i,idxRxns(j))) = 1;
                        d = d+1;
                    else
                        idxv = getIndexVariable(model,'v',i,idxRxns(j));
                        Aineq(t,idxv) = 1/(model.Kcat_f(idxRxns(j),l)*model.epsilon);
                    end
                end
            end
            t = t+1;                
        end
    end
            
    Aineq = [Aineq;Aineq2];
    bineq = [bineq,bineq2];
end
