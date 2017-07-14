function [Aeq,beq,Aineq,bineq,lb,ub] = formulateConstraints(model)
% formulateConstraints returns the equality and inequality constraint matrices and
% right hand side vectors, and the lower and upper bounds needed to run
% deFBA given a deFBA model structure
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
% Aeq           matrix containing all equality constraints for deFBA: rows = constraints, columns = LP variables
% beq           vector containing the right hand side of all equality constraints in Aeq
% Aineq         matrix containing all inequality constraints for deFBA: rows = constraints, columns = LP variables
% bineq         vector containing the right hand side of all inequality constraints in Aineq, s.t. Aineq*x<=bineq, where x is the solution vector
% lb            vector containing all variable lower bounds
% ub            vector containing all variable upper bounds
% 
% Alexandra Reimers 14/07/2017

    [Aeq1,beq1] = getConstraintStorageAndExtracellularDynamics(model);
    disp('formulated ydot=Sv') 
    
    [Aeq2,beq2] = getConstraintMacromoleculeDynamics(model);
    disp('formulated pdot=Sv')
    
    [Aeq3,beq3] = getConstraintMetaboliteSteadyState(model);
    disp('formulated constraint steady state')
    
    [AineqH,bineqH] = getConstraintEnzymeCapacity(model);
    disp('formulated constraint kcat')    
    
    [AeqH,beqH] = getConstraintDiscretization(model);
    disp('formulated constraint yi=yi-1+h/2ydot')
    
    [Aeq6,beq6] = getConstraintInitialAmounts(model);
    disp('formulated constraint initial y and bound on initial p')
    
    [Aineq6,bineq6] = getConstraintBiomassComposition(model);
    disp('formulated biomass composition constraint')
    
    [AineqM,bineqM] = getConstraintMaintenance(model);
    disp('formulated maintenance constraint')
    
    [lb,ub] = getBounds(model);
    disp('formulated variable bounds')
    
    Aeq = [Aeq1;Aeq2;Aeq3;AeqH;Aeq6];
    beq = [beq1,beq2,beq3,beqH,beq6];
    
    Aineq = [AineqH;Aineq6;AineqM];
    bineq  = [bineqH,bineq6,bineqM];
end