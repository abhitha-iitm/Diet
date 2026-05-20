%% =========================================
% OUTPUT FOLDER
%% =========================================

saveFolder = 'FVA_results_ROI';
if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end

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
% LOAD ROI
%% =========================================

roiTable = readtable('rois.xlsx');
roiRxns = cellstr(roiTable{:,1});

fprintf('Loaded %d ROI reactions\n', length(roiRxns));

%% =========================================
% 🔴 BASELINE FVA (NO DIET)
%% =========================================

fprintf('\n===== BASELINE FVA (NO DIET) =====\n');

model_CT = readCbModel('./Multitissue_Models/MulModel_CT.mat');
model_PC = readCbModel('./Multitissue_Models/MulModel_PC.mat');

model_CT = convert_EX_to_diet(model_CT);
model_PC = convert_EX_to_diet(model_PC);

models = {model_CT, model_PC};

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

model_CT = models{1};
model_PC = models{2};

% ROI filtering
roi_valid = intersect(roiRxns, model_CT.rxns);
roi_valid = intersect(roi_valid, model_PC.rxns);

fprintf('Baseline valid ROI: %d\n', length(roi_valid));

[min_CT_base, max_CT_base] = fluxVariability(model_CT, 100, 'max', roi_valid);
[min_PC_base, max_PC_base] = fluxVariability(model_PC, 100, 'max', roi_valid);

span_CT_base = max_CT_base - min_CT_base;
span_PC_base = max_PC_base - min_PC_base;

save(fullfile(saveFolder,'FVA_ROI_BASELINE.mat'), ...
    'roi_valid', ...
    'min_CT_base','max_CT_base','span_CT_base', ...
    'min_PC_base','max_PC_base','span_PC_base');

fprintf('✅ BASELINE DONE\n');

%% =========================================
% 🔵 LOOP OVER DIETS
%% =========================================

for d = 1:length(dietFiles)

    fprintf('\n===========================\n');
    fprintf('Diet: %s\n', dietFiles{d});
    fprintf('===========================\n');

    dietName = erase(dietFiles{d}, '.tsv');
    dietPath = fullfile('Diets', dietFiles{d});

    %% LOAD CLEAN MODEL

    model_CT = readCbModel('./Multitissue_Models/MulModel_CT.mat');
    model_PC = readCbModel('./Multitissue_Models/MulModel_PC.mat');

    model_CT = convert_EX_to_diet(model_CT);
    model_PC = convert_EX_to_diet(model_PC);

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

    %% LOAD OPTIMIZED MODEL

    data = load(fullfile('results-final-correct', ...
        ['nutritionResults_' dietName '.mat']));

    model_PC_opt = data.newDietModel;

    model_PC_opt.c(:) = 0;
    for i = 1:length(objRxns)
        idx = find(strcmp(model_PC_opt.rxns, objRxns{i}));
        if ~isempty(idx)
            model_PC_opt.c(idx) = weights(i);
        end
    end
    model_PC_opt.osenseStr = 'max';

    %% ROI FILTERING (CRITICAL)

    roi_valid = intersect(roiRxns, model_CT_diet.rxns);
    roi_valid = intersect(roi_valid, model_PC_diet.rxns);
    roi_valid = intersect(roi_valid, model_PC_opt.rxns);

    fprintf('Valid ROI reactions: %d\n', length(roi_valid));

    %% FVA (DIET)

    fprintf('Running FVA (diet)...\n');

    [min_CT, max_CT] = fluxVariability(model_CT_diet, 100, 'max', roi_valid);
    [min_PC, max_PC] = fluxVariability(model_PC_diet, 100, 'max', roi_valid);

    %% FVA (OPTIMIZED)

    fprintf('Running FVA (optimized)...\n');

    [min_PCopt, max_PCopt] = fluxVariability(model_PC_opt, 100, 'max', roi_valid);

    %% SPANS

    span_CT     = max_CT     - min_CT;
    span_PC     = max_PC     - min_PC;
    span_PCopt  = max_PCopt  - min_PCopt;

    %% SAVE

    save(fullfile(saveFolder, ['FVA_ROI_' dietName '.mat']), ...
        'roi_valid', ...
        'min_CT','max_CT','span_CT', ...
        'min_PC','max_PC','span_PC', ...
        'min_PCopt','max_PCopt','span_PCopt');

    fprintf('Completed: %s\n', dietName);

end

fprintf('\n✅ FULL ROI-FVA COMPLETED SUCCESSFULLY\n');