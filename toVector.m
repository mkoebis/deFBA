function w = toVector(str, model)
%toVector transform a struct of strtraint coefficients or solution struct into vector
%
% w = toVector(str,model)
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
% the order of variables is
% y0, v, ydot, pdot, y, p, p0, vRev
%
%INPUTS
% str        Struct that contains the coefficients for the constraints or the solution entries.
%            The struct should have the fields: y0, v, ydot, pdot, y, p, p0, vRev
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
% w           vector version of the coefficients or of the solution
%
% A note on scaling: For improved numerical stability, we scale some of the coefficients
%   v(internal reactions) - Unit: mmol/*h
%   v(exchange reactions) - Unit: mmol/*h
%   v(quota reactions) 	  - Unit: 1/*h
%   v(enzyme synthesis)   - Unit: mmol*epsilon/h
%   p(quota)              - Unit: pg
%   p(enzymes)            - Unit: epsilon * mmol
%   pdot(quota)           - Unit: pg/h
%   pdot(enzymes)         - Unit: epsilon * mmol/h
%   p0(quota)             - Unit: pg
%   p0(enzymes)           - Unit: epsilon * mmol
%   E                     - Unit: epsilon * mmol
%   y                     - Unit: mmol
%   ydot                  - Unit: mmol/h
%   y0                    - Unit: mmol
%
% Alexandra Reimers  17/07/2017
  
    % check structure sizes for consistency
    if (length(str.y0)~=model.sizeYmet || size(str.v,1)~=model.N || size(str.v,2)~=model.noRxn ||...
            size(str.ydot,1)~=model.N || size(str.ydot,2)~=model.sizeYmet || ...
            size(str.pdot,1)~=model.N || size(str.pdot,2)~=model.sizePmet || ...
            size(str.y,1)~=model.N || size(str.y,2)~=model.sizeYmet || ...
            size(str.p,1)~=model.N || size(str.p,2)~=model.sizePmet || ...
            length(str.p0)~=model.sizePmet)
        
        msgID = 'toVector:BadSize';
        msg = 'Structure sizes incompatible with the model';
        baseException = MException(msgID,msg);
        throw(baseException)
    end


    y0 = str.y0;
    v = reshape(str.v,1,numel(str.v));
    ydot = reshape(str.ydot,1,numel(str.ydot));
    pdot = reshape(str.pdot,1,numel(str.pdot));
    y = reshape(str.y,1,numel(str.y));
    p = reshape(str.p,1,numel(str.p));
    p0 = reshape(str.p0,1,numel(str.p0));
    vRev = reshape(str.vRev,1,numel(str.vRev));
    w = sparse([y0,v,ydot,pdot,y,p,p0,vRev]);    
end
