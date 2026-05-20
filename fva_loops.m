%% =========================================
% DIETS
%% =========================================

dietFiles = {
    %'EU.tsv'
    %'glutenFree.tsv'
    %'highprotein.tsv'
    %'type2diabetes.tsv'
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
% LOOP OVER DIETS
%% =========================================

for d = 1:length(dietFiles)

    fprintf('\n===========================\n');
    fprintf('Diet: %s\n', dietFiles{d});
    fprintf('===========================\n');

    dietName = erase(dietFiles{d}, '.tsv');
    dietPath = fullfile('Diets', dietFiles{d});

    %% =====================================
    % 🔴 LOAD CLEAN MODEL EACH TIME (CRITICAL FIX)
    %% =====================================

    model_CT = readCbModel('./Multitissue_Models/MulModel_CT.mat');
    model_PC = readCbModel('./Multitissue_Models/MulModel_PC.mat');

    model_CT = convert_EX_to_diet(model_CT);
    model_PC = convert_EX_to_diet(model_PC);

    %% =====================================
    % APPLY DIET
    %% =====================================

    [model_CT_diet, ~, ~] = setDietBoundsFromFile(model_CT, dietPath);
    [model_PC_diet, ~, ~] = setDietBoundsFromFile(model_PC, dietPath);

    %% =====================================
    % SET OBJECTIVE
    %% =====================================

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

    [min_CT_full, max_CT_full] = fluxVariability(model_CT_diet, 100, 'max', model_CT_diet.rxns);
    [min_PC_full, max_PC_full] = fluxVariability(model_PC_diet, 100, 'max', model_PC_diet.rxns);

    % ALIGN CT vs PC_diet
    [commonRxns, idx_CT, idx_PC] = intersect(model_CT_diet.rxns, model_PC_diet.rxns, 'stable');

    min_CT_aligned = min_CT_full(idx_CT);
    max_CT_aligned = max_CT_full(idx_CT);

    min_PC_aligned = min_PC_full(idx_PC);
    max_PC_aligned = max_PC_full(idx_PC);

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

    [min_PCopt_full, max_PCopt_full] = fluxVariability(model_PC_opt, 100, 'max', model_PC_opt.rxns);

    % ALIGN CT vs PC_opt
    [commonRxns_opt, idx_CT2, idx_PCopt] = intersect(model_CT_diet.rxns, model_PC_opt.rxns, 'stable');

    min_CT_opt = min_CT_full(idx_CT2);
    max_CT_opt = max_CT_full(idx_CT2);

    min_PCopt_aligned = min_PCopt_full(idx_PCopt);
    max_PCopt_aligned = max_PCopt_full(idx_PCopt);

    %% =====================================
    % SAVE
    %% =====================================

    save(fullfile(saveFolder, ['FVA_' dietName '.mat']), ...
        'commonRxns', ...
        'min_CT_aligned','max_CT_aligned', ...
        'min_PC_aligned','max_PC_aligned', ...
        'commonRxns_opt', ...
        'min_CT_opt','max_CT_opt', ...
        'min_PCopt_aligned','max_PCopt_aligned');

    fprintf('Completed: %s\n', dietName);

end

fprintf('\n✅ FINAL FVA COMPLETED SUCCESSFULLY\n');
%%
clear; clc;

dietNames = {'High Protein','Type2 Diabetes','Mediterranean','DACH','High Fiber',...
    'Vegan','Unhealthy','High Fat Low Carb','Gluten Free','Vegetarian','EU'};

before = [42.40,43.24,42.90,42.57,42.40,42.24,44.74,44.41,43.41,41.24,43.41];
after  = [36.39,37.80,37.73,37.70,37.85,37.90,40.57,40.90,40.40,39.40,42.07];

figure('Color','w','Position',[100 100 1100 480])

barData = [before; after]';
b = bar(barData,'grouped','BarWidth',0.75);

% Professional colors
b(1).FaceColor = [0.6 0.6 0.6];   % baseline
b(2).FaceColor = [0.2 0.45 0.75]; % optimized

% X-axis
xticks(1:length(dietNames))
xticklabels(dietNames)
xtickangle(35)

ylabel('Dysregulation (%)','FontSize',12)

% 🔴 FULL MEANINGFUL LEGEND
legend({
    'Baseline disease state (BASE vs PC)', ...
    'After diet optimization (BASE vs PCOPT)'
    },...
    'Location','northoutside',...
    'Orientation','horizontal',...
    'Box','off',...
    'FontSize',10)

% Axis styling
set(gca,'FontSize',11,'LineWidth',1.2)
box off
ylim([34 46])

title('Effect of Dietary Optimization on Metabolic Dysregulation',...
    'FontWeight','bold','FontSize',13)

% Adjust spacing (prevents overlap)
ax = gca;
ax.Position = [0.06 0.18 0.9 0.65];

exportgraphics(gcf,'Figure1_Final.png','Resolution',300)