%% =========================================
% INITIALIZE COBRA
%% =========================================

%initCobraToolbox(false);
%changeCobraSolver('gurobi','LP');

%% =========================================
% LOAD BASE MODELS
%% =========================================
if isempty(gcp('nocreate'))
    parpool('local');
end

changeCobraSolverParams('LP','Threads',8);

model_CT = readCbModel('./Multitissue_Models/MulModel_CT.mat');
model_PC = readCbModel('./Multitissue_Models/MulModel_PC.mat');

model_CT = convert_EX_to_diet(model_CT);
model_PC = convert_EX_to_diet(model_PC);

%% =========================================
% DIETS
%% =========================================

dietFiles = {
    'EU.tsv'
    'glutenFree.tsv'
    'highprotein.tsv'
    'type2diabetes.tsv'
    'mediterranean.tsv'
    'highfiber.tsv'
    'unhealthy.tsv'
    'vegan.tsv'
    'vegetarian.tsv'
    'DACH.tsv'
    'highfatlowcarb.tsv'
};

%% =========================================
% OUTPUT FOLDER
%% =========================================

saveFolder = 'FVA_results_final';
if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end

%% =========================================
% OBJECTIVE
%% =========================================

objRxns = {
    'SK_ATPtm'
    'SK_FAOXC204'
    'AD_ACCOAC'
    'SK_biomass_maintenance'
    'AD_biomass_maintenance'
    'GN_P450SCC1m'
    'GN_biomass_maintenance'
    'OO_ATPtm'
    'OO_biomass_maintenance'
    'EN_SERPT'
    'EN_SMS'
    'EN_biomass_maintenance'
};

weights = [
    0.23;0.09;0.26;0.01;0.01;
    0.13;0.01;0.11;0.01;0.06;0.07;0.01];

%% =========================================
% 🔴 BASELINE FVA (NO DIET)
%% =========================================

fprintf('\n===== BASELINE FVA =====\n');

model_CT_base = model_CT;
model_PC_base = model_PC;

models = {model_CT_base, model_PC_base};

for m = 1:2
    models{m}.c(:) = 0;

    for i = 1:length(objRxns)
        idx = find(strcmp(models{m}.rxns, objRxns{i}));
        if ~isempty(idx)
            models{m}.c(idx) = weights(i);
        end
    end

    models{m}.osenseStr = 'max';
end

model_CT_base = models{1};
model_PC_base = models{2};

[min_CT_base, max_CT_base] = fluxVariability(model_CT_base, 100, 'max', model_CT_base.rxns);
[min_PC_base, max_PC_base] = fluxVariability(model_PC_base, 100, 'max', model_PC_base.rxns);

[commonRxns_base, idx_CT_b, idx_PC_b] = intersect(model_CT_base.rxns, model_PC_base.rxns, 'stable');

min_CT_base = min_CT_base(idx_CT_b);
max_CT_base = max_CT_base(idx_CT_b);

min_PC_base = min_PC_base(idx_PC_b);
max_PC_base = max_PC_base(idx_PC_b);

save(fullfile(saveFolder,'FVA_baseline.mat'), ...
    'commonRxns_base', ...
    'min_CT_base','max_CT_base', ...
    'min_PC_base','max_PC_base');

%% =========================================
% LOOP OVER DIETS
%% =========================================

for d = 1:length(dietFiles)

    fprintf('\n===========================\n');
    fprintf('Diet: %s\n', dietFiles{d});
    fprintf('===========================\n');

    dietName = erase(dietFiles{d}, '.tsv');
    dietPath = fullfile('Diets', dietFiles{d});

    %% APPLY DIET

    [model_CT_diet, ~, ~] = setDietBoundsFromFile(model_CT, dietPath);
    [model_PC_diet, ~, ~] = setDietBoundsFromFile(model_PC, dietPath);

    %% SET OBJECTIVE

    models = {model_CT_diet, model_PC_diet};

    for m = 1:2
        models{m}.c(:) = 0;

        for i = 1:length(objRxns)
            idx = find(strcmp(models{m}.rxns, objRxns{i}));
            if ~isempty(idx)
                models{m}.c(idx) = weights(i);
            end
        end

        models{m}.osenseStr = 'max';
    end

    model_CT_diet = models{1};
    model_PC_diet = models{2};

    %% =====================================
    % 🔵 FVA (WITH DIET)
    %% =====================================

    fprintf('Running FVA (diet)...\n');

    [min_CT, max_CT] = fluxVariability(model_CT_diet, 100, 'max', model_CT_diet.rxns);
    [min_PC, max_PC] = fluxVariability(model_PC_diet, 100, 'max', model_PC_diet.rxns);

    [commonRxns, idx_CT, idx_PC] = intersect(model_CT_diet.rxns, model_PC_diet.rxns, 'stable');

    min_CT = min_CT(idx_CT);
    max_CT = max_CT(idx_CT);

    min_PC = min_PC(idx_PC);
    max_PC = max_PC(idx_PC);

    %% =====================================
    % 🔵 LOAD OPTIMIZED MODEL
    %% =====================================

    data = load(fullfile('results-final-correct', ...
        ['nutritionResults_' dietName '.mat']));

    model_PC_opt = data.newDietModel;

    % set objective for optimized model
    model_PC_opt.c(:) = 0;

    for i = 1:length(objRxns)
        idx = find(strcmp(model_PC_opt.rxns, objRxns{i}));
        if ~isempty(idx)
            model_PC_opt.c(idx) = weights(i);
        end
    end

    model_PC_opt.osenseStr = 'max';

    %% =====================================
    % 🔵 FVA (OPTIMIZED)
    %% =====================================

    fprintf('Running FVA (optimized)...\n');

    [min_PCopt, max_PCopt] = fluxVariability(model_PC_opt, 100, 'max', model_PC_opt.rxns);

    % ALIGN with CT_diet (IMPORTANT)
    [commonRxns_opt, idx_CT2, idx_PCopt] = intersect(model_CT_diet.rxns, model_PC_opt.rxns, 'stable');

    min_CT_opt = min_CT(idx_CT2);
    max_CT_opt = max_CT(idx_CT2);

    min_PCopt = min_PCopt(idx_PCopt);
    max_PCopt = max_PCopt(idx_PCopt);

    %% =====================================
    % SAVE
    %% =====================================

    save(fullfile(saveFolder, ['FVA_' dietName '.mat']), ...
        'commonRxns', ...
        'min_CT','max_CT','min_PC','max_PC', ...
        'commonRxns_opt', ...
        'min_CT_opt','max_CT_opt', ...
        'min_PCopt','max_PCopt');

end

fprintf('\n✅ FINAL FVA COMPLETED SUCCESSFULLY\n');