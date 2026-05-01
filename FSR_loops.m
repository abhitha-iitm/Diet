%% =========================================
% DIET LIST
%% =========================================

dietNames = {
    'EU'
    'glutenFree'
    'highprotein'
    'type2diabetes'
    'mediterranean'
    'highfiber'
    'unhealthy'
    'vegan'
    'vegetarian'
    'DACH'
    'highfatlowcarb'
};

%% =========================================
% PATHS
%% =========================================

inputFolder  = 'FVA_results_final';
outputFolder = 'FSR_results_final';

if ~exist(outputFolder,'dir')
    mkdir(outputFolder)
end

%% =========================================
% LOOP OVER DIETS
%% =========================================

for d = 1:length(dietNames)

    diet = dietNames{d};
    fprintf('\n===========================\n');
    fprintf('Processing diet: %s\n', diet);
    fprintf('===========================\n');

    filePath = fullfile(inputFolder, ['FVA_' diet '.mat']);

    if ~isfile(filePath)
        fprintf('Skipping %s (file not found)\n', diet);
        continue;
    end

    data = load(filePath);

    %% =====================================
    % CT vs PC (DIET)
    %% =====================================

    rxns = string(data.commonRxns(:));

    min_CT = data.min_CT;
    max_CT = data.max_CT;

    min_PC = data.min_PC;
    max_PC = data.max_PC;

    %% =====================================
    % CT vs PC_opt (OPTIMIZED)
    %% =====================================

    rxns_opt = string(data.commonRxns_opt(:));

    min_CT_opt = data.min_CT_opt;
    max_CT_opt = data.max_CT_opt;

    min_PCopt = data.min_PCopt;
    max_PCopt = data.max_PCopt;

    %% =====================================
    % RUN BOTH COMPARISONS
    %% =====================================

    comparisons = {
        rxns,     min_CT,     max_CT,     min_PC,     max_PC,     ['CT_vs_PC_' diet];
        rxns_opt, min_CT_opt, max_CT_opt, min_PCopt,  max_PCopt,  ['CT_vs_PCOPT_' diet];
    };

    for c = 1:size(comparisons,1)

        rxnList = comparisons{c,1};
        min_normal = comparisons{c,2};
        max_normal = comparisons{c,3};
        min_disease = comparisons{c,4};
        max_disease = comparisons{c,5};
        tag = comparisons{c,6};

        fprintf('Running FSR: %s\n', tag)

        %% =====================================
        % FSR CALCULATION (WITH SAFE TOL)
        %% =====================================

        CT_span = max_normal - min_normal;
        Disease_span = max_disease - min_disease;

        nRxns = length(rxnList);
        FSR = NaN(nRxns,1);

        tol = 1e-9;

        for i = 1:nRxns
            denom = CT_span(i);

            if abs(denom) > tol
                FSR(i) = Disease_span(i) / denom;
            else
                FSR(i) = NaN;   % safer than 0
            end
        end

        %% =====================================
        % CLASSIFICATION
        %% =====================================

        Status = repmat("Unchanged", nRxns, 1);

        Status(FSR >= 2) = "Upregulated";
        Status(FSR <= 0.5) = "Downregulated";

        %% =====================================
        % STATS
        %% =====================================

        num_up = sum(FSR >= 2, 'omitnan');
        num_down = sum(FSR <= 0.5, 'omitnan');

        fprintf('Total: %d | Up: %d | Down: %d\n', ...
            nRxns, num_up, num_down);

        %% =====================================
        % TABLE
        %% =====================================

        FSR_table = table(rxnList, FSR, Status, ...
            'VariableNames', {'Reaction','FSR','RegulationStatus'});

        %% =====================================
        % SAVE
        %% =====================================

        writetable(FSR_table, ...
            fullfile(outputFolder, ['FSR_' tag '.csv']));

        writetable(FSR_table(Status~="Unchanged",:), ...
            fullfile(outputFolder, ['FSR_Dysregulated_' tag '.csv']));
    end
end

fprintf('\n✅ FSR COMPLETED SUCCESSFULLY\n');