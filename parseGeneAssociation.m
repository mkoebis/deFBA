% parse the model.rules to figure out how each reaction is catalyzed
% INPUT    
% model     model structure that has a field model.rules
%
% OUTPUT   
% ands      array of size length(model.rxns) that tells for each reaction
%           how many '&' are in the rules entry
% ors       array of size length(model.rxns) that tells for each reaction
%           how many '|' are in the rules entry
% We assume for simplicity that each rules has either only and signs or
% only or signs and throw a warning if that is not the case


function [ands,ors,geneCombs] = parseGeneAssociation(model)
    ands = zeros(length(model.rxns),1);
    ors = zeros(length(model.rxns),1);
    geneCombs = {};
    for i=1:length(model.rules)
        % find ands
        found = strfind(model.rules(i),'&');
        ands(i) = length(found{1});
        % find ors
        found = strfind(model.rules(i),'|');
        ors(i) = length(found{1});
        
        % remove parentheses and add to geneCombs
        if ors(i)>0
            % count number of unique conjunctions
            mysplit = strsplit(model.grRules{i},' OR ');
            for j=1:length(mysplit)
                mysplit{j} = strrep(mysplit{j},'(','');
                mysplit{j} = strrep(mysplit{j},')','');
                if ~strcmp(mysplit{j},'')
                    geneCombs{length(geneCombs)+1} = mysplit{j};
                end
            end
        else
            mysplit = model.grRules{i};
            mysplit = strrep(mysplit,'(','');
            mysplit = strrep(mysplit,')','');
            if ~strcmp(mysplit,'')
                geneCombs{length(geneCombs)+1} = mysplit;
            end
        end
     end
end