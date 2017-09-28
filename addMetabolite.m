function [ newmodel ] = addMetabolite(model,metID,metName)
%addMetabolite adds a Metabolite to the model
%
% newModel = addMetabolite(model,metID,metName)
%
% metID and metName string arguments either have to be a single value or Cell
% Arrays. Charge and b have to be double arrays
%
%INPUTS
% model         Cobra model structure
% metID         The ID(s) of the metabolite(s) (will be the identifier in model.mets)
%
%OPTIONAL INPUTS
% metName       Human readable name(s) (String)
%
%OUTPUT
% newModel      COBRA model with added metabolite(s)
%
% Alexandra Reimers 26/09/2017

varName={'char','cell','numeric','logical'};

if ~isempty(metID)
    if ~isa(metID,'cell')
        for i=1:numel(varName);
            if isa(metID,varName{i});
                type=i;
            end
            
        end
        if type==1;
            metID = {metID};
        else
            %errorMsg=sprintf('The type of the metID should be ''char'',but here the provided metID is ''%d''',varName{type});
            errorMsg=varName{type};
            
            errorMsg=['The type of the provided metID should be ''char'' or ''cell'', but here the provided metID is ', errorMsg]; 
            errordlg(errorMsg);            
        end
    end
else
    errordlg('metID is empty');
end

if nargin < 3
    metName = cell(1,numel(metID));
    metName(:) = {''};
else
    if ~isa(metName,'cell')
        if ~isa(metName,'char')
            fprintf('Wrong Argument class for metName: %s ; should be char or cell\n',class(metName));
            return
        else
            metName = {metName};
        end
    end
    if numel(metID) ~= numel(metName)
        fprintf('Inconsistent Argument length metID (%i) and metName(%i)\n',numel(metID),numel(metName));
        return
    end
end

for i = 1:numel(metID)
    cmetID = metID{i};
    if isempty(find(ismember(model.mets,cmetID)))
        model.S(end+1,:) = 0;
        model.mets{end+1} = cmetID;
        if (isfield(model,'metNames'))      %Prompts to add missing info if desired
            cmetName = metName{i};
            if strcmp(cmetName,'')
                model.metNames{end+1,1} = regexprep(cmetID,'(\[.+\]) | (\(.+\))','') ;                
                warning(['Metabolite name for ' metID{i} ' set to ' model.metNames{end}]);
            else
                model.metNames{end+1,1} = metName{i} ;                
            end
        
        end
    end
end

newmodel = model;

end

