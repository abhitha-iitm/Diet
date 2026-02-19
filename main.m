%initCobraToolbox(false);
% changeCobraSolver('gurobi','LP');

%% LOAD MODELS
model_CT = readCbModel('./Multitissue_Models/MulModel_CT.mat');
model_PC = readCbModel('./Multitissue_Models/MulModel_PC.mat');

model_CT = convert_EX_to_diet(model_CT);
model_PC = convert_EX_to_diet(model_PC);

%% LIST OF DIETS TO TEST
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
};

resultsFolder = 'results-sameDiets';
if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

%% LOAD ROI LIST
roiTable = readtable('rois.xlsx');
roi_list = roiTable{:,1};

%% LOOP OVER DIETS
for d = 1:length(dietFiles)

    fprintf('Running diet: %s\n', dietFiles{d});

    dietPath = fullfile('Diets', dietFiles{d});

    %% APPLY SAME DIET TO BOTH MODELS
    [model_CT_diet, ~, ~] = setDietBoundsFromFile(model_CT, dietPath);
    [model_PC_diet, ~, ~] = setDietBoundsFromFile(model_PC, dietPath);

    %% SET WEIGHTED BIOMASS OBJECTIVE
    tissueBiomassRxns = {
        'SK_biomass_maintenance'
        'AD_biomass_maintenance'
        'GN_biomass_maintenance'
        'OO_biomass_maintenance'
        'EN_biomass_maintenance'
    };

    weights = [0.2 0.2 0.2 0.2 0.2];

    model_CT_diet.c(:) = 0;
    model_PC_diet.c(:) = 0;

    for i = 1:length(tissueBiomassRxns)
        idx1 = find(strcmp(model_CT_diet.rxns,tissueBiomassRxns{i}));
        idx2 = find(strcmp(model_PC_diet.rxns,tissueBiomassRxns{i}));

        if ~isempty(idx1)
            model_CT_diet.c(idx1) = weights(i);
        end
        if ~isempty(idx2)
            model_PC_diet.c(idx2) = weights(i);
        end
    end

    %% OPTIMIZE CONTROL MODEL
    sol_CT = optimizeCbModel(model_CT_diet,'max');
    fprintf('CT objective: %.4f\n', sol_CT.f);

    %% OPTIMIZE PCOS MODEL
    sol_PC = optimizeCbModel(model_PC_diet,'max');
    fprintf('PCOS objective: %.4f\n', sol_PC.f);

    %% EXTRACT ROI FLUXES
    control_flux = zeros(length(roi_list),1);
    pcos_flux    = zeros(length(roi_list),1);

    for i = 1:length(roi_list)
        idxCT = strcmp(model_CT_diet.rxns, roi_list{i});
        idxPC = strcmp(model_PC_diet.rxns, roi_list{i});

        control_flux(i) = sol_CT.v(idxCT);
        pcos_flux(i)    = sol_PC.v(idxPC);
    end

    %% SET OPTIONS FOR NUTRITION ALGORITHM
    options.pcosFlux   = pcos_flux;
    options.targetFlux = control_flux;
    options.roiWeights = ones(1,length(roi_list));
    options.display    = 'on';

    % allow only foods present in diet
    allowedDietRxns = model_CT_diet.rxns( ...
        model_CT_diet.lb < 0 & startsWith(model_CT_diet.rxns,'Diet_EX_'));

    options.targetedDietRxns = ...
        [allowedDietRxns, num2cell(ones(length(allowedDietRxns),1))];

    %% RUN NUTRITION ALGORITHM
    [newDietModel, pointsModel, roiFlux, pointsModelSln, ...
        menuChanges] = ...
        nutritionAlgorithm_new(model_PC_diet, roi_list, options);

    %% DISPLAY RESULTS
    fprintf('\nSuggested dietary changes:\n');
    disp(menuChanges);

    %% SAVE RESULTS
    dietName = erase(dietFiles{d},'.tsv');

    save(fullfile(resultsFolder, ...
        ['nutritionResults_' dietName '.mat']), ...
        'newDietModel', 'pointsModel', 'roiFlux', ...
        'pointsModelSln', 'menuChanges', ...
        'control_flux', 'pcos_flux');

    %% SAVE MENU CHANGES AS CSV
    writetable(menuChanges, ...
        fullfile(resultsFolder, ...
        ['menuChanges_' dietName '.csv']));

end

fprintf('\nAll diets completed.\n');
