
%initCobraToolbox(false);
%changeCobraSolver('gurobi','LP');

%% LOADING MODELS
model_CT = readCbModel('./Multitissue_Models/MulModel_CT.mat');
model_PC = readCbModel('./Multitissue_Models/MulModel_PC.mat');

model_CT = convert_EX_to_diet(model_CT);
model_PC = convert_EX_to_diet(model_PC);

%% LIST OF DIETS
dietFiles = {
    'mediterranean.tsv'
    };
%% =========================================
% DEFINE ALL WEIGHT CASES
%% =========================================
tissueBiomassRxns = {
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
allWeights = {
    % CASE 1: Equal (0.2 each)
    % [
    %     0.095;0.095;0.195;0.01;0.01;
    %     0.19;0.01;0.19;0.01;0.095;0.095;0.01
    % ],

    % CASE 2: 0.5 / 0.5
    % [
    %     0.18;0.08;0.20;0.01;0.01;
    %     0.17;0.01;0.15;0.01;0.09;0.08;0.01
    % ],

    % CASE 3: 0.6 / 0.4
    [
        0.23;0.09;0.26;0.01;0.01;
        0.13;0.01;0.11;0.01;0.06;0.07;0.01
    ],

    % CASE 4: 0.7 / 0.3
    % [
    %     0.26;0.1;0.28;0.01;0.01;
    %     0.11;0.01;0.09;0.01;0.05;0.06;0.01
    % ],

    % CASE 5: PCOS
    % [
    %     0.2;0.08;0.18;0.01;0.01;
    %     0.17;0.01;0.16;0.01;0.08;0.08;0.01
    % ]
};

caseNames = {
    'case1_equal'
    'case2_50_50'
    'case3_60_40'
    'case4_70_30'
    'case5_pcos'
};

%% LOADING ROI LIST
roiTable = readtable('rois.xlsx');
roi_list = roiTable{:,1};

%% =========================================
% LOOP OVER CASES
%% =========================================

for c = 1:length(allWeights)

    weights = allWeights{c};
    caseFolder = fullfile('', caseNames{c});

    if ~exist(caseFolder,'dir')
        mkdir(caseFolder);
    end

    fprintf('\n===========================\n');
    fprintf('Running %s\n', caseNames{c});
    fprintf('===========================\n');

    %% LOOP OVER DIETS
    for d = 1:length(dietFiles)

        fprintf('Running diet: %s\n', dietFiles{d});

        dietPath = fullfile('Diets', dietFiles{d});

        %% APPLY DIET
        [model_CT_diet, ~, ~] = setDietBoundsFromFile(model_CT, dietPath);
        [model_PC_diet, ~, ~] = setDietBoundsFromFile(model_PC, dietPath);

        %% RESET OBJECTIVE
        model_CT_diet.c(:) = 0;
        model_PC_diet.c(:) = 0;

        %% APPLY WEIGHTS
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

        %% OPTIMIZATION
        sol_CT = optimizeCbModel(model_CT_diet,'max');
        sol_PC = optimizeCbModel(model_PC_diet,'max');

        %% ROI FLUXES
        control_flux = zeros(length(roi_list),1);
        pcos_flux    = zeros(length(roi_list),1);

        for i = 1:length(roi_list)
            idxCT = strcmp(model_CT_diet.rxns, roi_list{i});
            idxPC = strcmp(model_PC_diet.rxns, roi_list{i});

            control_flux(i) = sol_CT.v(idxCT);
            pcos_flux(i)    = sol_PC.v(idxPC);
        end

        %% OPTIONS
        options.pcosFlux   = pcos_flux;
        options.targetFlux = control_flux;
        options.roiWeights = ones(1,length(roi_list));
        options.display    = 'on';

        allowedDietRxns = model_CT_diet.rxns( ...
            model_CT_diet.lb < 0 & startsWith(model_CT_diet.rxns,'Diet_EX_'));

        options.targetedDietRxns = ...
            [allowedDietRxns, num2cell(ones(length(allowedDietRxns),1))];

        %% NUTRITION ALGO
        [newDietModel, pointsModel, roiFlux, pointsModelSln, menuChanges] = ...
            nutritionAlgorithm_new(model_PC_diet, roi_list, options);

        %% SAVE
        dietName = erase(dietFiles{d},'.tsv');

        save(fullfile(caseFolder, ...
            ['nutritionResults_' dietName '.mat']), ...
            'newDietModel','pointsModel','roiFlux', ...
            'pointsModelSln','menuChanges', ...
            'control_flux','pcos_flux');

        writetable(menuChanges, ...
            fullfile(caseFolder, ...
            ['menuChanges_' dietName '.csv']));
    end
end