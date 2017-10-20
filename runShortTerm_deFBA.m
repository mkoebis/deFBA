% runShortTerm_deFBA runs short term deFBA as described in
% https://arxiv.org/abs/1609.08961
% on the model structure. 
% The prediction horizon is the one set in the field deFBAmodel.tf
% Please make sure that the order of the reactions and metabolites in the 
% stoichiometric matrix follows the pattern:
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
%   tf                          final time of simulation a.k.a. prediction horizon (hours)
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
% noIterations    The number of iteration steps. If you e.g. want to run
%                 the short term deFBA for a period of 1 hour, then make sure that
%                 noIteration*deFBAmodel.tf/deFBAmodel.N = 1h.
% solver          String containing the name of the solver that should be
%                 used. Current options are 'cplex' or 'soplex'. For cplex to work
%                 the cplex matlab directory should be in the matlab path. For soplex to
%                 work, the solver needs to be installed and the used needs to edit the
%                 path to the soplex executable in the initSolver.m file.
%
%OUTPUT
% result          solution structure. Note that it may not be possible to
%                 use toVector and turn this structure into a solution vector compatible
%                 with deFBA because of number of time steps missmatch
%
% Alexandra Reimers 02/10/2017

function result = runShortTerm_deFBA(deFBAmodel,noIterations,solver)
    % first iteration separate, to set up and get p0
    sol0 = run_deFBA(deFBAmodel,solver);
    strsol = toStruct(sol0,deFBAmodel);

    % save result for this iteration
    result.p = [strsol.p0;strsol.p(1,:)];
    result.y = [strsol.y0;strsol.y(1,:)];
    result.v = strsol.v(1,:);
    result.pdot = strsol.pdot(1,:);
    result.ydot = strsol.ydot(1,:);    
    
    % set starting point for next iteration
    deFBAmodel.Y0 = result.y(end,:);
    deFBAmodel.p0 = result.p(end,:);
    
    % then regular iterations till max number of iterations reached
    for i = 1:noIterations-1
        % solve regular deFBA
        sol = run_deFBA(deFBAmodel,solver);
        strsol = toStruct(sol,deFBAmodel);
        
        % save result for this iteration
        result.p = [result.p;strsol.p(1,:)];
        result.y = [result.y;strsol.y(1,:)];
        result.v = [result.v;strsol.v(1,:)];
        result.pdot = [result.pdot;strsol.pdot(1,:)];
        result.ydot = [result.ydot;strsol.ydot(1,:)];    

        % set starting point for next iteration
        deFBAmodel.Y0 = result.y(end,:);
        deFBAmodel.p0 = result.p(end,:);
        
        % move switches one step closer and if they are at zero implement
        % and remove them
        toRemove = [];
        if  isfield(deFBAmodel,'switches')
            for j=1:length(deFBAmodel.switches)
                deFBAmodel.switches(j).n = deFBAmodel.switches(j).n-1;

                if deFBAmodel.switches(j).n==0
                    deFBAmodel.Y0(deFBAmodel.switches(j).idx) = deFBAmodel.switches(j).amount;
                    toRemove = [toRemove,j];
                end
            end
            deFBAmodel.switches = deFBAmodel.switches(setdiff(1:1:length(deFBAmodel.switches),toRemove));
        end
        
        % move flux bounds (but not those of knockouts) 
        if isfield(deFBAmodel,'rxnExtraBounds')
            if ~isempty(deFBAmodel.rxnExtraBounds)
                if isfield(deFBAmodel.rxnExtraBounds{1},'rxnID')
                    for j=1:length(deFBAmodel.rxnExtraBounds)
                        if deFBAmodel.rxnExtraBounds{j}.KOflag == false
                            toRemove = [];
                            deFBAmodel.rxnExtraBounds{j}.timePoint = deFBAmodel.rxnExtraBounds{j}.timePoint-1;

                            if any(deFBAmodel.rxnExtraBounds{j}.timePoint==0)
                                toRemove = find(deFBAmodel.rxnExtraBounds{j}.timePoint==0);
                            end
                            deFBAmodel.rxnExtraBounds{j}.timePoint = deFBAmodel.rxnExtraBounds{j}.timePoint(setdiff(1:1:length(deFBAmodel.rxnExtraBounds{j}.timePoint),toRemove));
                        end
                    end
                    % remove empty bounds
                    idx = false(length(deFBAmodel.rxnExtraBounds),1);
                    for k=1:length(deFBAmodel.rxnExtraBounds)
                        idx(k) = isempty(deFBAmodel.rxnExtraBounds{k}.timePoint);
                    end
                    
                    deFBAmodel.rxnExtraBounds = deFBAmodel.rxnExtraBounds(~idx);
                    if isempty(deFBAmodel.rxnExtraBounds)
                        deFBAmodel = rmfield(deFBAmodel,'rxnExtraBounds');
                    end
                end
            end
        end
        
        if norm(deFBAmodel.Y0(deFBAmodel.noStorage+1:end))<1e-5
            break;
        end
    end
end
