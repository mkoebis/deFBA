function outmodel = writedeFBAModel(model,fileName,compSymbolList,compNameList,externalCompID)
%writedeFBAModel Write a deFBA model to SBML using the ram standard
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
% fileName          File name for output file
% compSymbolList    List of compartment ids
% compNameList      List of compartment names corresponding to compSymbolList
% externalCompID    ID of the external compartment
%
% OPTIONAL OUTPUTS
% outmodel          The sbml structure.
% 
% Alexandra Reimers 02/08/2017

    outmodel = writeSBML_ram(model,fileName,compSymbolList,compNameList,externalCompID);
    
    % remove double gene products (libsbml bug)
    fid1 = fopen(fileName,'r'); 
    filename2 = sprintf('%s_new.xml',fileName(1:end-4));
    fid2 = fopen(filename2,'w'); 
    while ~feof(fid1)
        line = fgets(fid1); %# read line by line
        if isempty(strfind(line,'fbc:id="gp_'))
            fprintf(fid2,'%s',line); %# write the line to the new file
        end
    end
    fclose(fid1);
    fclose(fid2);
        
    % correct infix problem with gene associations (unnecessary extra gene products, libsbml bug)
    sbmlText = fileread(filename2);
    sbmlText = strrep(sbmlText,'gp_','');
    sbmlText = strrep(sbmlText,'__ZERO__','0');
    sbmlText = strrep(sbmlText,'__ONE__','1');
    sbmlText = strrep(sbmlText,'__TWO__','2');
    sbmlText = strrep(sbmlText,'__THREE__','3');
    sbmlText = strrep(sbmlText,'__FOUR__','4');
    sbmlText = strrep(sbmlText,'__FIVE__','5');
    sbmlText = strrep(sbmlText,'__SIX__','6');
    sbmlText = strrep(sbmlText,'__SEVEN__','7');
    sbmlText = strrep(sbmlText,'__EIGHT__','8');
    sbmlText = strrep(sbmlText,'__NINE__','9');
    fid = fopen(fileName,'w');
    fprintf(fid,'%s',sbmlText);
    fclose(fid);
    s = sprintf('%s_new.xml',fileName(1:end-4));
    delete(s);
end
