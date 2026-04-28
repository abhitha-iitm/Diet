% %% =========================================
% % INITIALIZE COBRA
% %% =========================================

% initCobraToolbox(false)
% changeCobraSolver('gurobi','LP')

% %% =========================================
% % DIET LIST
% %% =========================================

% dietNames = {
%     'EU'
%     'glutenFree'
%     'highprotein'
%     'type2diabetes'
%     'mediterranean'
%     'highfiber'
%     'unhealthy'
%     'vegan'
%     'vegetarian'
%     };

% %% =========================================
% % OBJECTIVE DEFINITIONS
% %% =========================================

% tissueAndBiomassRxns = { ...
%     'SK_ATPtm', ...
%     'SK_biomass_maintenance', ...
%     'AD_ACCOAC', ...
%     'AD_biomass_maintenance', ...
%     'GN_biomass_maintenance', ...
%     'GN_P450SCC1m', ...
%     'OO_ATPtm', ...
%     'OO_biomass_maintenance', ...
%     'EN_biomass_maintenance', ...
%     'EN_SERPT', ...
%     'EN_SMS' ...
%     };

% tissueAndBiomassWeights=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.0666,0.0667,0.0667];

% %% =========================================
% % CREATE FVA RESULTS FOLDER
% %% =========================================

% baseFolder = 'FVA_results';

% if ~exist(baseFolder,'dir')
%     mkdir(baseFolder)
% end

% %% =========================================
% % LOOP OVER CASES
% %% =========================================

% for caseID = 5

%     fprintf('\n===========================\n')
%     fprintf('Running Case %d\n',caseID)
%     fprintf('===========================\n')

%     caseFolder = 'Case5_biomassAndTissue_specific_sameDietControl';
%     objRxns = tissueAndBiomassRxns;
%     weights = tissueAndBiomassWeights;

%     casePath = fullfile(baseFolder,caseFolder);

%     if ~exist(casePath,'dir')
%         mkdir(casePath)
%     end

%     %% =====================================
%     % DIET LOOP
%     %% =====================================

%     for d = 1:length(dietNames)

%         diet = dietNames{d};

%         saveFile = fullfile(casePath, ['FVA_' diet '.mat']);
%         if exist(saveFile, 'file')
%             fprintf('⏭️ Skipping %s (already completed)\n', diet);
%             continue;
%         end

%         fprintf('\nProcessing diet: %s\n',diet)

%         %% LOAD MODELS

%         load(fullfile('Diet_models','control',['CT_' diet '.mat']))
%         model_control = model_CT_diet;

%         load(fullfile('Diet_models','pcos',['PC_' diet '.mat']))
%         model_PC = model_PC_diet;

%         data = load(fullfile('results-tso',['nutritionResults_' diet '.mat']));
%         model_PC_opt = data.newDietModel;

%         %% SET OBJECTIVES

%         models = {model_control, model_PC, model_PC_opt};

%         for m = 1:3
%             models{m}.c(:) = 0;
%             for i = 1:length(objRxns)
%                 idx = strcmp(models{m}.rxns,objRxns{i});
%                 models{m}.c(idx) = weights(i);
%             end
%             models{m}.osenseStr = 'max';
%         end

%         model_control = models{1};
%         model_PC = models{2};
%         model_PC_opt = models{3};

%         %% RUN FVA

%         fprintf('Running FVA: Control + %s\n',diet)
%         [min_control,max_control] = fluxVariability(model_control,100);

%         fprintf('Running FVA: PCOS + %s\n',diet)
%         [min_PC,max_PC] = fluxVariability(model_PC,100);

%         fprintf('Running FVA: PCOS + optimized %s\n',diet)
%         [min_PC_opt,max_PC_opt] = fluxVariability(model_PC_opt,100);

%         %% =====================================
%         % 🔥 FINAL ALIGNMENT (REFERENCE-BASED)
%         %% =====================================

%         % Reference (control order)
%         refRxns = regexprep(lower(strtrim(cellstr(model_control.rxns))), '[^a-z0-9_]', '');

%         pcRxns  = regexprep(lower(strtrim(cellstr(model_PC.rxns))),     '[^a-z0-9_]', '');
%         optRxns = regexprep(lower(strtrim(cellstr(model_PC_opt.rxns))), '[^a-z0-9_]', '');

%         [tfP, locP] = ismember(refRxns, pcRxns);
%         [tfO, locO] = ismember(refRxns, optRxns);

%         min_PC_aligned      = nan(size(refRxns));
%         max_PC_aligned      = nan(size(refRxns));
%         min_PC_opt_aligned  = nan(size(refRxns));
%         max_PC_opt_aligned  = nan(size(refRxns));

%         % fill
%         min_PC_aligned(tfP)      = min_PC(locP(tfP));
%         max_PC_aligned(tfP)      = max_PC(locP(tfP));

%         min_PC_opt_aligned(tfO)  = min_PC_opt(locO(tfO));
%         max_PC_opt_aligned(tfO)  = max_PC_opt(locO(tfO));

%         % final assignment
%         rxns = refRxns;
%         min_PC = min_PC_aligned;
%         max_PC = max_PC_aligned;
%         min_PC_opt = min_PC_opt_aligned;
%         max_PC_opt = max_PC_opt_aligned;

%         assert(numel(min_control) == numel(refRxns), 'Control FVA length mismatch');
% assert(numel(min_PC) == numel(refRxns), 'PC FVA length mismatch');
% assert(numel(min_PC_opt) == numel(refRxns), 'PC opt FVA length mismatch');


%         %% SAVE

%         save(saveFile,...
%             'rxns',...
%             'min_control','max_control',...
%             'min_PC','max_PC',...
%             'min_PC_opt','max_PC_opt');

%         fprintf('Saved → %s\n', diet);

%     end
% end

% fprintf('\nFVA completed successfully.\n');

%% =========================================
% INITIALIZE COBRA
%% =========================================

%initCobraToolbox(false)
%changeCobraSolver('gurobi','LP')

if isempty(gcp('nocreate'))
    parpool('local');
end

%% =========================================
% DIET LIST
%% =========================================

dietNames = {
    %'EU'
    % 'glutenFree'
    % 'highprotein'
    % 'type2diabetes'
     'mediterranean'
    %'highfiber'
    %'unhealthy'
    % 'vegan'
    % 'vegetarian'
    };

%% =========================================
% NEW OBJECTIVE DEFINITIONS
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
allWeights = {
    % CASE 1: Equal (0.2 each)
    % [
    %     0.095;0.095;0.195;0.01;0.01;
    %     0.19;0.01;0.19;0.01;0.095;0.095;0.01
    % ],

    % CASE 2: 0.5 / 0.5
    [
        0.18;0.08;0.20;0.01;0.01;
        0.17;0.01;0.15;0.01;0.09;0.08;0.01
    ],

    % CASE 3: 0.6 / 0.4
    [
        0.23;0.09;0.26;0.01;0.01;
        0.13;0.01;0.11;0.01;0.06;0.07;0.01
    ],

    % % CASE 4: 0.7 / 0.3
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
    % 'case1_equal'
    'case2_50_50'
    'case3_60_40'
    % 'case4_70_30'
    % 'case5_pcos'
};

%% =========================================
% CREATE FOLDER
%% =========================================

baseFolder = 'FVA_results-M-allCases';

if ~exist(baseFolder,'dir')
    mkdir(baseFolder)
end

%% =========================================
% LOOP OVER CASES
%% =========================================

for caseID = 1:length(allWeights)
    weights = allWeights{caseID};
    fprintf('\n===========================\n')
    fprintf('Running Case %d\n',caseID)
    fprintf('===========================\n')

    caseFolder = caseNames{caseID};
    casePath = fullfile(baseFolder,caseFolder);

    if ~exist(casePath,'dir')
        mkdir(casePath)
    end

    %% =====================================
    % DIET LOOP
    %% =====================================

    for d = 1:length(dietNames)

        diet = dietNames{d};

        saveFile = fullfile(casePath, ['FVA_' diet '.mat']);
        if exist(saveFile, 'file')
            fprintf('Skipping %s (already completed)\n', diet);
            continue;
        end

        fprintf('\nProcessing diet: %s\n',diet)

        %% LOAD MODELS

        load(fullfile('Diet_models','control',['CT_' diet '.mat']))
        model_control = model_CT_diet;

        load(fullfile('Diet_models','pcos',['PC_' diet '.mat']))
        model_PC = model_PC_diet;

        data = load(fullfile('NA-results_diffWeights_M', caseNames{caseID}, ...
    ['nutritionResults_' diet '.mat']));
        model_PC_opt = data.newDietModel;

        %% SET OBJECTIVES (SAFE VERSION)

        models = {model_control, model_PC, model_PC_opt};

        for m = 1:3
            models{m}.c(:) = 0;

            for i = 1:length(objRxns)
                idx = find(strcmp(models{m}.rxns, objRxns{i}));

                if ~isempty(idx)
                    models{m}.c(idx) = weights(i);
                end
            end

            models{m}.osenseStr = 'max';
        end

        model_control = models{1};
        model_PC      = models{2};
        model_PC_opt  = models{3};

        %% =====================================
        % RUN FAST FVA
        %% =====================================

        fprintf('Running FVA (CT)...\n');
        [min_CT, max_CT] = fluxVariability(model_control, 100, 'max', model_control.rxns);

        fprintf('Running FVA (PC)...\n');
        [min_PC, max_PC] = fluxVariability(model_PC, 100, 'max', model_PC.rxns);

        fprintf('Running FVA (PC_opt)...\n');
        [min_PCopt, max_PCopt] = fluxVariability(model_PC_opt, 100, 'max', model_PC_opt.rxns);

        %% =====================================
        % COMMON REACTIONS (CT vs PC)
        %% =====================================

        [commonRxns, idx_CT, idx_PC] = intersect(model_control.rxns, model_PC.rxns, 'stable');

        min_CT_common = min_CT(idx_CT);
        max_CT_common = max_CT(idx_CT);

        min_PC_common = min_PC(idx_PC);
        max_PC_common = max_PC(idx_PC);

        %% =====================================
        % COMMON REACTIONS (CT vs PC_opt)
        %% =====================================

        [commonRxns_opt, idx_CT2, idx_PCopt] = intersect(model_control.rxns, model_PC_opt.rxns, 'stable');

        min_CT_common_opt = min_CT(idx_CT2);
        max_CT_common_opt = max_CT(idx_CT2);

        min_PCopt_common = min_PCopt(idx_PCopt);
        max_PCopt_common = max_PCopt(idx_PCopt);

        %% =====================================
        % SAVE RESULTS
        %% =====================================

save(saveFile, ...
    'commonRxns', ...
    'min_CT_common', 'max_CT_common', ...
    'min_PC_common', 'max_PC_common', ...
    'commonRxns_opt', ...
    'min_CT_common_opt', 'max_CT_common_opt', ...
    'min_PCopt_common', 'max_PCopt_common');

    end
end

