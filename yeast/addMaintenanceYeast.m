function model = addMaintenanceYeast(model)
    model = addReaction(model,'ATP_maintenance',{'s_0434[c03]','s_0803[c03]','s_0394[c03]','s_1322[c03]','s_0794[c03]'},[-1 -1 1 1 1],false,'');
end