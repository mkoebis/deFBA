% incorporate into deFBA model kcat corresponding to rxns
% assumption: forward and reverse direction have the same kcat for
% reversible reactions
% for reactions not in rxns assign median of kcat
function deFBAmodel = assignDataKcatYeast(deFBAmodel,rxns,kcat)
    % assignment for reactions with known kcat
    indicesRest = true(deFBAmodel.noRxn,1);
    for i=1:length(rxns)
        idx = find(strcmp(deFBAmodel.rxnECNumbers,rxns{i}));
        % we ignore kcat values that are too small, because they can be wrong
        if ~isempty(idx) && kcat(i)>0.1
            for k=1:length(idx)
                indicesRest(idx(k))=false;
                enzidx = find(deFBAmodel.Kcat_f(idx(k),:));
                for j=1:length(enzidx)
                    deFBAmodel.Kcat_f(idx(k),enzidx(j)) = kcat(i)*3600;
                end
            end
        end
    end
    % assignment for the rest
    m = 10*3600;% the moderately efficient enzyme
    for i=1:length(indicesRest)
        if indicesRest(i)
            idx = i;
            enzidx = find(deFBAmodel.Kcat_f(idx,:));
            for j=1:length(enzidx)
                deFBAmodel.Kcat_f(idx,enzidx(j)) = m;
            end
        end
    end
end