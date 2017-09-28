function formulas = printRxnFormula(model,rxnAbbrList,metNameFlag)
%formulas = printRxnFormula(model,rxnAbbrList,printFlag,lineChangeFlag,metNameFlag,fid,directionFlag)
%printRxnFormula Print the reaction formulas for a list of reactions
%
%
%INPUTS
% model             COBRA model structure
%
%OPTIONAL INPUTS
% rxnAbbrList       Abbrs of reactions whose formulas are to be printed
% metNameFlag       print full met names instead of abbreviations 
%                   (Default = false)
%
%OUTPUT
% formulas          Cell array containing formulas of specified reactions
%
% Alexandra Reimers 27/09/2017

if (nargin < 2)
    rxnAbbrList = model.rxns;
end


lineChangeFlag = true;

if (nargin <3)
    metNameFlag = true;
end

fid = 1;


if (~iscell(rxnAbbrList))
    if (strcmp(rxnAbbrList,'all'))
        rxnAbbrList = model.rxns;
    else
        rxnAbbrTmp = rxnAbbrList;
        clear rxnAbbrList;
        rxnAbbrList{1} = rxnAbbrTmp;
    end
end

for i = 1:length(rxnAbbrList);

    rxnAbbr = rxnAbbrList{i};

    rxnID = findRxnIDs(model,rxnAbbr);

    fprintf(fid,'%s\t',rxnAbbr);
    
    if (rxnID > 0)

        Srxn = full(model.S(:,rxnID));

        Sprod = (Srxn(Srxn > 0));
        if metNameFlag
            prodMets = model.metNames(Srxn > 0);
        else
            prodMets = model.mets(Srxn > 0);
        end
        
        Sreact = (Srxn(Srxn < 0));
        if metNameFlag
            reactMets = model.metNames(Srxn < 0);
        else
            reactMets = model.mets(Srxn < 0);
        end
        
        formulaStr = '';
        for j = 1:length(reactMets)
            if (j > 1)
                fprintf(fid,'+ ');
                formulaStr = [formulaStr '+ '];
            end
            if (abs(Sreact(j)) ~= 1)
                fprintf(fid,'%f %s ',abs(Sreact(j)),reactMets{j});
                formulaStr = [formulaStr num2str(abs(Sreact(j))) ' ' reactMets{j} ' '];
            else
                fprintf(fid,'%s ',reactMets{j});                
                formulaStr = [formulaStr reactMets{j} ' '];
            end
        end

        if (model.rev(rxnID))
            fprintf(fid,'\t<=>\t');
            formulaStr = [formulaStr ' <=> '];
        else
            fprintf(fid,'\t->\t');            
            formulaStr = [formulaStr ' -> '];
        end
        
       
        for j = 1:length(prodMets)
            if (j > 1)
                fprintf(fid,'+ ');
                formulaStr = [formulaStr '+ '];
            end
            if (Sprod(j) ~= 1)
                fprintf(fid,'%f %s ',Sprod(j),prodMets{j});
                formulaStr = [formulaStr num2str(Sprod(j)) ' ' prodMets{j} ' '];
            else
                fprintf(fid,'%s ',prodMets{j});                
                formulaStr = [formulaStr prodMets{j} ' '];
            end
        end
        
    else
        fprintf(fid,'not in model');
        formulaStr = 'NA';
    end
   
    if (lineChangeFlag)
        fprintf(fid,'\n');
    end    
    formulas{i} = formulaStr;
end
formulas = formulas';
