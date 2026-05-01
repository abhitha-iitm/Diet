% %initCobraToolbox(false);
% %changeCobraSolver('gurobi','LP');

% % LOADING MODELS
% model_CT = readCbModel('./Multitissue_Models/MulModel_CT.mat');
% model_PC = readCbModel('./Multitissue_Models/MulModel_PC.mat');

% model_CT = convert_EX_to_diet(model_CT);
% model_PC = convert_EX_to_diet(model_PC);

% % LIST OF DIETS
% dietFiles = {
%     'EU.tsv'
%     'glutenFree.tsv'
%     'highprotein.tsv'
%     'type2diabetes.tsv'
%     'mediterranean.tsv'
%     'highfiber.tsv'
%     'unhealthy.tsv'
%     'vegan.tsv'
%     'vegetarian.tsv'
%     'DACH.tsv'
%     'highfatlowcarb.tsv'
%     };

% % [model_CT_diet, ~, ~] = setDietBoundsFromFile(model_CT, fullfile('Diets',dietFiles{1}));
% resultsFolder = 'results-new-objectives-newlimits';
% if ~exist(resultsFolder,'dir')
%     mkdir(resultsFolder);
% end

% % 🔴 ADD THIS (store scores)
% numDiets = length(dietFiles);
% dysregulationScores = zeros(numDiets,1);

% % LOADING ROI LIST
% roiTable = readtable('rois.xlsx');
% roi_list = roiTable{:,1};

% % LOOP OVER DIETS
% for d = 1:length(dietFiles)

%     fprintf('Running diet: %s\n', dietFiles{d});

%     dietPath = fullfile('Diets', dietFiles{d});

%     % APPLYING SAME DIET TO BOTH MODELS
%     [model_CT_diet, ~, ~] = setDietBoundsFromFile(model_CT, dietPath);
%     [model_PC_diet, ~, ~] = setDietBoundsFromFile(model_PC, dietPath);

%     % SET WEIGHTED BIOMASS OBJECTIVE
%     % tissueBiomassRxns = {
%     %     'SK_biomass_maintenance'
%     %     'AD_biomass_maintenance'
%     %     'GN_biomass_maintenance'
%     %     'OO_biomass_maintenance'
%     %     'EN_biomass_maintenance'
%     %     };
%     % weights = [0.2 0.2 0.2 0.2 0.2];
%     % tissueBiomassRxns = {'SK_ATPtm','AD_ACCOAC',...
%     %                  'GN_biomass_maintenance','GN_P450SCC1m',...
%     %                  'OO_ATPtm','EN_SERPT','EN_SMS'};

%     % weights=[0.2,0.2,0.1,0.1,0.2,0.1,0.1];

%     % tissueBiomassRxns = { ...
%     %     'SK_ATPtm', ...
%     %     'SK_biomass_maintenance', ...
%     %     'AD_ACCOAC', ...
%     %     'AD_biomass_maintenance', ...
%     %     'GN_biomass_maintenance', ...
%     %     'GN_P450SCC1m', ...
%     %     'OO_ATPtm', ...
%     %     'OO_biomass_maintenance', ...
%     %     'EN_biomass_maintenance', ...
%     %     'EN_SERPT', ...
%     %     'EN_SMS' ...
%     %     };

%     % weights=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.0666,0.0667,0.0667];

%     tissueBiomassRxns = {
%         'SK_ATPtm'
%         'SK_FAOXC204'
%         'AD_ACCOAC'
%         'SK_biomass_maintenance'
%         'AD_biomass_maintenance'

%         'GN_P450SCC1m'
%         'GN_biomass_maintenance'
%         'OO_ATPtm'
%         'OO_biomass_maintenance'
%         'EN_SERPT'
%         'EN_SMS'
%         'EN_biomass_maintenance'
%         };

%     weights =     [
%         0.23;0.09;0.26;0.01;0.01;
%         0.13;0.01;0.11;0.01;0.06;0.07;0.01
%         ];

%     model_CT_diet.c(:) = 0;
%     model_PC_diet.c(:) = 0;

%     for i = 1:length(tissueBiomassRxns)
%         idx1 = find(strcmp(model_CT_diet.rxns,tissueBiomassRxns{i}));
%         idx2 = find(strcmp(model_PC_diet.rxns,tissueBiomassRxns{i}));

%         if ~isempty(idx1)
%             model_CT_diet.c(idx1) = weights(i);
%         end
%         if ~isempty(idx2)
%             model_PC_diet.c(idx2) = weights(i);
%         end
%     end

%     % OPTIMIZE CONTROL MODEL
%     sol_CT = optimizeCbModel(model_CT_diet,'max');
%     fprintf('CT objective: %.4f\n', sol_CT.f);

%     % OPTIMIZE PCOS MODEL
%     sol_PC = optimizeCbModel(model_PC_diet,'max');
%     fprintf('PCOS objective: %.4f\n', sol_PC.f);

%     % EXTRACT ROI FLUXES
%     control_flux = zeros(length(roi_list),1);
%     pcos_flux    = zeros(length(roi_list),1);

%     for i = 1:length(roi_list)
%         idxCT = strcmp(model_CT_diet.rxns, roi_list{i});
%         idxPC = strcmp(model_PC_diet.rxns, roi_list{i});

%         control_flux(i) = sol_CT.v(idxCT);
%         pcos_flux(i)    = sol_PC.v(idxPC);
%     end

%     %% SET OPTIONS FOR NUTRITION ALGORITHM
%     options.pcosFlux   = pcos_flux;
%     options.targetFlux = control_flux;
%     options.roiWeights = ones(1,length(roi_list));
%     options.display    = 'on';

%     % allow only foods present in diet
%     allowedDietRxns = model_CT_diet.rxns( ...
%         model_CT_diet.lb < 0 & startsWith(model_CT_diet.rxns,'Diet_EX_'));

%     options.targetedDietRxns = ...
%         [allowedDietRxns, num2cell(ones(length(allowedDietRxns),1))];

%     % CORRECT DEBUG BLOCK
%     dietIdx = find(contains(model_PC_diet.rxns,'Diet_EX'));

%     lbVals = model_PC_diet.lb(dietIdx);

%     activeIdx = dietIdx(lbVals < 0);

%     fprintf('\n===== DIET SUMMARY (%s) =====\n', dietFiles{d})

%     % Count of active metabolites
%     fprintf('Number of active diet metabolites: %d\n', length(activeIdx))

%     % Total uptake
%     fprintf('Total uptake: %f\n', sum(lbVals))

%     % Show active ones
%     disp(table(model_PC_diet.rxns(activeIdx), ...
%         model_PC_diet.lb(activeIdx)))
%     %% 🔴 BASELINE PCOS deviation

% dev_before = sum(abs(pcos_flux - control_flux));

% fprintf('\n🟡 Baseline Deviation (%s): %f\n', dietFiles{d}, dev_before);

%     %% NUTRITION ALGORITHM
%     [newDietModel, pointsModel, roiFlux, pointsModelSln, ...
%         menuChanges] = ...
%         nutritionAlgorithm_new(model_PC_diet, roi_list, options);

%     %%  DYSREGULATION SCORE (ROI-based — BEST FOR YOU)

%     score = pointsModelSln.f;
%     dysregulationScores(d) = score;

% fprintf('\n🔴 Dysregulation Score (%s): %f\n', dietFiles{d}, score);
%     %% RESULTS
%     fprintf('\nSuggested dietary changes:\n');
%     disp(menuChanges);

%     %% RESULTS
%     dietName = erase(dietFiles{d},'.tsv');

%     save(fullfile(resultsFolder, ...
%         ['nutritionResults_' dietName '.mat']), ...
%         'newDietModel', 'pointsModel', 'roiFlux', ...
%         'pointsModelSln', 'menuChanges', ...
%         'control_flux', 'pcos_flux');

%     %% MENU CHANGES AS CSV
%     writetable(menuChanges, ...
%         fullfile(resultsFolder, ...
%         ['menuChanges_' dietName '.csv']));

% end

% % FINAL COMPARISON TABLE
% fprintf('\n===== FINAL DIET COMPARISON =====\n');

% T = table(dietFiles, dysregulationScores);
% disp(T);

% %  SORT (BEST → WORST)
% [sortedScores, idx] = sort(dysregulationScores);

% fprintf('\n===== SORTED (BEST → WORST) =====\n');
% disp(table(dietFiles(idx), sortedScores));

% % BAR PLOT (IMPORTANT FOR THESIS)
% figure;
% bar(dysregulationScores);
% xticks(1:numDiets);
% xticklabels(dietFiles);
% xtickangle(45);
% ylabel('Dysregulation Score');
% title('Diet Comparison (PCOS vs Control)');






%initCobraToolbox(false);
%changeCobraSolver('gurobi','LP');

% LOADING MODELS
model_CT = readCbModel('./Multitissue_Models/MulModel_CT.mat');
model_PC = readCbModel('./Multitissue_Models/MulModel_PC.mat');

model_CT = convert_EX_to_diet(model_CT);
model_PC = convert_EX_to_diet(model_PC);

% LIST OF DIETS
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

resultsFolder = 'results-final-correct';
if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

numDiets = length(dietFiles);

% 🔴 STORE ALL THREE SCORES
baselineScores   = zeros(numDiets,1);
afterDietScores  = zeros(numDiets,1);
finalScores      = zeros(numDiets,1);

% ROI LIST
roiTable = readtable('rois.xlsx');
roi_list = roiTable{:,1};

%% =========================
% 🔴 TRUE BASELINE (NO DIET)
%% =========================

fprintf('\n===== TRUE BASELINE (NO DIET) =====\n');

sol_CT_base = optimizeCbModel(model_CT,'max');
sol_PC_base = optimizeCbModel(model_PC,'max');

control_flux_base = zeros(length(roi_list),1);
pcos_flux_base    = zeros(length(roi_list),1);

for i = 1:length(roi_list)
    control_flux_base(i) = sol_CT_base.v(strcmp(model_CT.rxns, roi_list{i}));
    pcos_flux_base(i)    = sol_PC_base.v(strcmp(model_PC.rxns, roi_list{i}));
end

baseline_true = sum(abs(pcos_flux_base - control_flux_base));

fprintf('Baseline (no diet): %f\n', baseline_true);

%% =========================
% LOOP OVER DIETS
%% =========================

for d = 1:length(dietFiles)

    fprintf('\n===========================\n');
    fprintf('Running diet: %s\n', dietFiles{d});
    fprintf('===========================\n');

    dietPath = fullfile('Diets', dietFiles{d});

    [model_CT_diet, ~, ~] = setDietBoundsFromFile(model_CT, dietPath);
    [model_PC_diet, ~, ~] = setDietBoundsFromFile(model_PC, dietPath);

    %% SET OBJECTIVE
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

    weights = [
        0.23;0.09;0.26;0.01;0.01;
        0.13;0.01;0.11;0.01;0.06;0.07;0.01];

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

    %% SOLVE MODELS (AFTER DIET)
    sol_CT = optimizeCbModel(model_CT_diet,'max');
    sol_PC = optimizeCbModel(model_PC_diet,'max');

    control_flux = zeros(length(roi_list),1);
    pcos_flux    = zeros(length(roi_list),1);

    for i = 1:length(roi_list)
        control_flux(i) = sol_CT.v(strcmp(model_CT_diet.rxns, roi_list{i}));
        pcos_flux(i)    = sol_PC.v(strcmp(model_PC_diet.rxns, roi_list{i}));
    end

    %% 🟢 AFTER DIET SCORE (CORRECT WAY)
    score_diet = sum(abs(pcos_flux - control_flux));
    afterDietScores(d) = score_diet;

    fprintf('🟢 After Diet Score: %f\n', score_diet);

    %% OPTIONS FOR NUTRITION ALGORITHM
    options.pcosFlux   = pcos_flux;
    options.targetFlux = control_flux;
    options.roiWeights = ones(1,length(roi_list));
    options.display    = 'on';

    %% RUN NUTRITION ALGORITHM
    [newDietModel, pointsModel, roiFlux, pointsModelSln, ...
        menuChanges] = ...
        nutritionAlgorithm_new(model_PC_diet, roi_list, options);

    %% 🔴 FINAL SCORE (OPTIMIZED)
    score_final = pointsModelSln.f;
    finalScores(d) = score_final;

    fprintf('🔴 After Optimization Score: %f\n', score_final);

    %% SAVE
    dietName = erase(dietFiles{d},'.tsv');

    save(fullfile(resultsFolder, ...
        ['nutritionResults_' dietName '.mat']), ...
        'newDietModel','pointsModel','roiFlux', ...
        'pointsModelSln','menuChanges', ...
        'control_flux','pcos_flux');

    writetable(menuChanges, ...
        fullfile(resultsFolder, ...
        ['menuChanges_' dietName '.csv']));

end

%% =========================
% FINAL COMPARISON
%% =========================

fprintf('\n===== FINAL COMPARISON =====\n');

T = table(dietFiles, ...
          repmat(baseline_true,numDiets,1), ...
          afterDietScores, ...
          finalScores, ...
          'VariableNames',{'Diet','Baseline','AfterDiet','AfterOptimization'});

disp(T);

%% SORT
[sortedScores, idx] = sort(finalScores);

fprintf('\n===== SORTED (BEST → WORST) =====\n');
disp(table(dietFiles(idx), sortedScores));

%% PLOT
figure;
bar(finalScores);
xticklabels(dietFiles);
xtickangle(45);
ylabel('Dysregulation Score');
title('Optimized Diet Comparison');


%%
% Data (use your variables directly if already in workspace)
dietNames = dietFiles;

afterDiet = afterDietScores;
afterOpt  = finalScores;

% Create figure
figure;

% Bar plot
barData = [afterDiet, afterOpt];
b = bar(barData);

% Labels
xticks(1:length(dietNames));
xticklabels(dietNames);
xtickangle(45);

ylabel('Dysregulation Score');
title('Diet vs Optimized Diet Comparison');

% Legend
legend({'After Diet','After Optimization'}, 'Location','northwest');

% Optional: improve appearance
grid on;