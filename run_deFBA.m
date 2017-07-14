function result = run_deFBA(model,solver)
% run_deFBA runs deFBA on the model structure. Please make sure that the
% order of the reactions and metabolites in the stoichiometric matrix
% follows the pattern:
% 
% xxsyqqeer
% xxsyqqeer
% sssyqqeer
% yyyyqqeer
% qqqqqqeer
% qqqqqqeer
% eeeeeeeer
% eeeeeeeer
% rrrrrrrrr
% 
% where x corresponds to internal mets and rxns, 
%       s corresponds to storage and its production/consumption reactions,
%       y corresponds to extracellular nutrients and their production/consumption reactions,
%       q corresponds to quota compounds and their production reactions,
%       e corresponds to metabolic enzymes and their production reactions,
%       r corresponds to ribosome and its production reaction
%
%INPUT
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
% solver          String containing the name of the solver that should be
%                 used. Current options are 'cplex' or 'soplex'. For cplex to work
%                 the cplex matlab directory should be in the matlab path. For soplex to
%                 work, the solver needs to be installed and the used needs to edit the
%                 path to the soplex executable in the initSolver.m file.
%
%OUTPUT
% result          solution vector
%
% Alexandra Reimers 14/07/2017


    % model.noRxn might be supplied by the model, but we don't want to rely on it, because it is easily computable
    model.noRxn = length(model.rxns);
    % sizeP might be supplied by the model, but we don't want to rely on it, because it is easily computable
    model.sizePmet = length(model.enz)+model.sizeQuotaMet;

    % check if the kcat matrices have the correct size
    assert(size(model.Kcat_f,1) == model.noRxn);
    assert(size(model.Kcat_f,2) == model.sizePmet-model.sizeQuotaMet);

    % check if the stoichiometric matrix has the correct size
    assert(size(model.S,1) == length(model.mets));
    assert(size(model.S,2) == length(model.rxns));
    
    % compute Xmet
    model.sizeXmet = size(model.S,1)-model.sizeYmet-model.sizePmet;
    model.Xmet = model.mets(1:model.sizeXmet);
    % it is horrible if we do some indexing mistakes, so make the user check if this is correct
    fprintf('first internal metabolite: %s\n', model.Xmet{1})
    fprintf('last internal metabolite: %s\n', model.Xmet{end})
    
    % verify enzymes
    for i=1:model.sizePmet-model.sizeQuotaMet
        assert(strcmp(model.mets{model.sizeXmet+model.sizeYmet+model.sizeQuotaMet+i}, model.enz{i}));
    end
    
    % compute Xrxn
    model.sizeXrxn = size(model.S,2)-model.sizeYrxn-model.sizePrxn;
    model.Xrxn = model.rxns(1:model.sizeXrxn);
    % it is horrible if we do some indexing mistakes, so make the user check if this is correct
    fprintf('first internal reaction: %s\n', model.Xrxn{1})
    fprintf('last internal reaction: %s\n', model.Xrxn{end})
    
    % compute Yrxn
    model.Yrxn = model.rxns(model.sizeXrxn+1:model.sizeXrxn+model.sizeYrxn);
    % it is horrible if we do some indexing mistakes, so make the user check if this is correct
    fprintf('first exchange reaction: %s\n', model.Yrxn{1})
    fprintf('last exchange reaction: %s\n', model.Yrxn{end})
    
    % compyte Prxn
    model.Prxn = model.rxns(model.sizeXrxn+model.sizeYrxn+1:end);
    assert(length(model.Prxn) == model.sizePrxn)
    fprintf('first enzyme synthesis reaction: %s\n', model.Prxn{1})
    fprintf('last enzyme synthesis reaction: %s\n', model.Prxn{end})
    
    [Aeq,beq,Aineq,bineq,lb,ub] = formulateConstraints(model);
    c = getObjectiveCoefficients(model);
    
    % create cplex object
    lpProb = Cplex;
    % set parameters
    lpProb.Param.simplex.tolerances.optimality.Cur = 1e-9;
    lpProb.Param.simplex.tolerances.feasibility.Cur = 1e-9;
    lpProb.Param.feasopt.tolerance.Cur = 1e-9;
    lpProb.Param.sifting.display.Cur = 0;
    lpProb.Param.emphasis.numerical.Cur = 1;
    lpProb.Param.simplex.display.Cur = 0;
    lpProb.Param.tune.display.Cur = 3;
    lpProb.Param.barrier.convergetol.Cur = 1e-12;
    lpProb.Param.barrier.display.Cur = 0;
    lpProb.Param.threads.Cur = 12;
    lpProb.Param.workmem.Cur = 2048;     
    lpProb.Param.read.scale.Cur = 0;
     
    % add bounds and objective
    lpProb.addCols(c',[],lb',ub');
    lpProb.Model.sense = 'maximize';
    
    lpProb.addRows([beq,-inf*ones(1,size(Aineq,1))]',[Aeq;Aineq],[beq,bineq]');
    
    if strcmp(solver,'cplex')
            % solve
            lpProb.solve();
            
            % get solution if it exists
            if lpProb.Solution.status == 1
                result = lpProb.Solution.x;
            end
            
            % if a dilution term was used for better numerics we scale back
            % the solution. Please note that if the dilution term is too
            % large the problem may be infeasible and reducing it might
            % help make it feasible again
            if model.beta > 0
                [result, ~] = transformDilutionResult(result,model);
            end
            [result,~] = transformEpsilonScaling(result,model);
            
    else if strcmp(solver,'soplex')
            ts = octave_random(0,intmax,1,1);
            f = sprintf('deFBALP%d.lp',ts);
            lpProb.writeModel(f);
            
            % solve with soplex with iterative refinement 
            soplexPath = initSolver;            
            syscall = sprintf('%s -f1e-13 -o1e-13 -v5 -x --readmode=1 --solvemode=2 -c deFBALP%d.lp>sol%d',soplexPath,ts,ts);            
            [s, ~] = system(syscall,'-echo');
            
            % delete lp file
            delete(f);
            
            if s ~= 0
                error('error calling soplex');
            end
            
            % read solution
            f = sprintf('sol%d',ts);
            if exist(f,'file')
                result = readSoPlexSolution(model,f);
            end
            
            % if a dilution term was used for better numerics we scale back
            % the solution. Please note that if the dilution term is too
            % large the problem may be infeasible and reducing it might
            % help make it feasible again
            if model.beta > 0
                [result, ~] = transformDilutionResult(result,model);
            end
            [result,~] = transformEpsilonScaling(result,model);            
            
        else
            error('Solver option unknown')
        end
    end
end