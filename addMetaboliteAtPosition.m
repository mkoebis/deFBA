function model = addMetaboliteAtPosition(model,metID,metName,position)
    if position>1
        model.S = [model.S(1:position-1,:);zeros(1,size(model.S,2));model.S(position:end,:)];
        model.mets = [model.mets(1:position-1);metID;model.mets(position:end)];
        model.metNames = [model.metNames(1:position-1);metName;model.metNames(position:end)];
    elseif position==1
        model.S = [zeros(1,size(model.S,2));model.S];
        model.mets = [;metID;model.mets];
        model.metNames = [metName;model.metNames];
    end
end