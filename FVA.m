
%% Load models
% Model_CT_n = load('models/new_models/MulModel_CT.mat'); 
% model_CT_EU = Model_CT_n.Model_CT_n;

% Model_PC_n = load('models/new_models/MulModel_PC.mat'); 
% model_PC_HF = Model_PC_n.Model_PC_n;

% Model_newDiet = load('models/new_models/newDietModel.mat'); 
% newDietModel = Model_newDiet.newDietModel;

%% Define tissue biomass reactions and weights
tissueBiomassRxns = {
    'SK_biomass_maintenance'
    'AD_biomass_maintenance'
    'GN_biomass_maintenance'
    'OO_biomass_maintenance'
    'EN_biomass_maintenance'
};

weights = [0.2, 0.2, 0.2, 0.2, 0.2];

%% Function to set weighted biomass objective
function model = setWeightedBiomassObjective(model, tissueBiomassRxns, weights)
    model.c(:) = 0;
    for i = 1:length(tissueBiomassRxns)
        rxnIdx = find(strcmp(model.rxns, tissueBiomassRxns{i}));
        if ~isempty(rxnIdx)
            model.c(rxnIdx) = weights(i);
        else
            fprintf('Warning: %s not found in model.rxns\n', tissueBiomassRxns{i});
        end
    end
end
objIdx = find(newDietModel.c ~= 0);

disp('Objective reaction indices:');
disp(objIdx);

disp('Objective reaction names:');
disp(newDietModel.rxns(objIdx));

disp('Objective coefficients:');
disp(newDietModel.c(objIdx));


%% Apply weighted biomass objective to all models
model_CT_EU = setWeightedBiomassObjective(model_CT_EU, tissueBiomassRxns, weights);
model_PC_HF = setWeightedBiomassObjective(model_PC_HF, tissueBiomassRxns, weights);
newDietModel = setWeightedBiomassObjective(optimizedModel, tissueBiomassRxns, weights);

%% Optimize each model to confirm feasibility
sol_CT = optimizeCbModel(model_CT_EU, 'max');
fprintf('CT + EU weighted biomass = %.6f\n', sol_CT.f);

sol_PC = optimizeCbModel(model_PC_HF, 'max');
fprintf('PC + HF weighted biomass = %.6f\n', sol_PC.f);

sol_new = optimizeCbModel(newDietModel, 'max');
fprintf('PC + Optimized HF weighted biomass = %.6f\n', sol_new.f);

%% ================================
% Run FVA
%% ================================

fprintf('\nRunning FVA on CT + EU model...\n');
[minFlux_CT, maxFlux_CT] = fluxVariability(model_CT_EU, 100, 'max', model_CT_EU.rxns);

fprintf('Running FVA on PC + HF model...\n');
[minFlux_PC, maxFlux_PC] = fluxVariability(model_PC_HF, 100, 'max', model_PC_HF.rxns);

fprintf('Running FVA on PC + Optimized HF model...\n');
[minFlux_new, maxFlux_new] = fluxVariability(newDietModel, 100, 'max', newDietModel.rxns);

%% ================================
% Save FVA results
%% ================================

fva_CT_table = table(model_CT_EU.rxns, minFlux_CT, maxFlux_CT, ...
    'VariableNames', {'Reaction', 'MinFlux', 'MaxFlux'});
writetable(fva_CT_table, 'FVA_CT_EU.csv');

fva_PC_table = table(model_PC_HF.rxns, minFlux_PC, maxFlux_PC, ...
    'VariableNames', {'Reaction', 'MinFlux', 'MaxFlux'});
writetable(fva_PC_table, 'FVA_PC_HF.csv');

fva_new_table = table(newDietModel.rxns, minFlux_new, maxFlux_new, ...
    'VariableNames', {'Reaction', 'MinFlux', 'MaxFlux'});
writetable(fva_new_table, 'FVA_PC_OptimizedHF.csv');

fprintf('\nFVA completed successfully. Files saved.\n');
