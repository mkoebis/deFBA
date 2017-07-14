function idx = getIndexVariable(model,field,i,k)
% getIndexVariable computes index for a given variable in the solution 
% vector given its name, timepoint (optional), and entry
%     
% solutions have the following variables
% y0   -  dimensions model.sizeYmet
% v    -  dimensions model.N x model.noRxn
% ydot -  dimensions model.N x model.sizeYmet
% pdot -  dimensions model.N x model.sizePmet
% y    -  dimensions model.N x model.sizeYmet
% p    -  dimensions model.N x model.sizePmet
% p0   -  dimensions model.sizePmet
% vRev -  dimensions model.N x sum(model.rev)
%
% For all variables slicing is supported
%
% The order of variables is
% y0, v, ydot, pdot, y, p, p0, vRev
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
% field             String containing the type of variable (options: 'y0', 'v', 'ydot', 'pdot', 'y', 'p', 'p0', 'vRev')
% i                 First dimension index (time point for dynamic variables or entry for y0 and p0)
% j                 Second dimension index. Only for dynamic variables, the entry of that variable
% 
%OUTPUT
% idx               The variable index in the solution.
%
% Alexandra Reimers  14/07/2017
    
    % we have initial amount variables for each external and storage metabolite amount    
    % we have a flux variable for each time step, for each reaction
    sizeV = model.N*model.noRxn;
    % we have a derivative of the external and storage metabolite amount
    % for every time step, and external metabolite or storage metabolite
    sizeYdot = model.N*model.sizeYmet;
    % we have a derivative for every macromolecule amount for every time step, and every macromolecule
    sizePdot = model.N*model.sizePmet;
    % we have an external metabolite amount for each time step and external metabolite
    sizeY = model.N*model.sizeYmet;
    % we have a macromolecule amount at each time step for each macromolecule
    sizeP = model.N*model.sizePmet;
    % we have initial macromolecule amount variables for time 0 for each macromolecule
    sizeP0 = model.sizePmet;
    
    % for each flux variable of a reversible reaction we have in addition a helper variable vRev that helps us get the absolute value

    % the order of variables is
    % y0, v, ydot, pdot, y, p, p0, vRev
    if strcmp(field,'y0')
        idx = i;
    end
    
    start = model.sizeYmet;
    
    if strcmp(field,'v')
        [I,K] = ndgrid(i,k);
        I = reshape(I,1,numel(I));
        K = reshape(K,1,numel(K));
        idx = start + sub2ind([model.N,model.noRxn],I,K);
    end
    
    if strcmp(field,'ydot')
        [I,K] = ndgrid(i,k);
        I = reshape(I,1,numel(I));
        K = reshape(K,1,numel(K));
        idx = start + sizeV + sub2ind([model.N,model.sizeYmet],I,K);
    end
    
    if strcmp(field,'pdot')
        [I,K] = ndgrid(i,k);
        I = reshape(I,1,numel(I));
        K = reshape(K,1,numel(K));
        idx = start + sizeV + sizeYdot + sub2ind([model.N,model.sizePmet],I,K);
    end
    
    if strcmp(field,'y')
        [I,K] = ndgrid(i,k);
        I = reshape(I,1,numel(I));
        K = reshape(K,1,numel(K));
        idx = start + sizeV + sizeYdot + sizePdot + sub2ind([model.N,model.sizeYmet],I,K);
    end
    
    if strcmp(field,'p')
        [I,K] = ndgrid(i,k);
        I = reshape(I,1,numel(I));
        K = reshape(K,1,numel(K));
        idx = start + sizeV + sizeYdot + sizePdot + sizeY + sub2ind([model.N,model.sizePmet],I,K);
    end 
    
    if strcmp(field,'p0')
        idx = start + sizeV + sizeYdot + sizePdot + sizeY + sizeP + i ;
    end
      
    if strcmp(field,'vRev')
        [I,K] = ndgrid(i,k);
        I = reshape(I,1,numel(I));
        K = reshape(K,1,numel(K));
        idx = start + sizeV + sizeYdot + sizePdot + sizeY + sizeP + sizeP0 + sub2ind([model.N,sum(model.rev)],I,K);
    end
end
