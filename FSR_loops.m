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
% % INPUT + OUTPUT BASE FOLDERS
% %% =========================================

% inputBase = 'FVA_results';
% outputBase = 'FSR_results';

% if ~exist(outputBase,'dir')
%     mkdir(outputBase)
% end

% %% =========================================
% % LOOP OVER CASES
% %% =========================================

% for caseID = 1:5

%     fprintf('\n===========================\n')
%     fprintf('Running FSR Case %d\n',caseID)
%     fprintf('===========================\n')

%     switch caseID
%         case 1
%             caseFolder = 'Case1_biomass_EUcontrol';
%         case 2
%             caseFolder = 'Case2_tissue_EUcontrol';
%         case 3
%             caseFolder = 'Case3_biomass_sameDietControl';
%         case 4
%             caseFolder = 'Case4_tissue_sameDietControl';
%         case 5
%             caseFolder = 'Case5_biomassAndTissue_specific_sameDietControl';
%     end

%     inputPath = fullfile(inputBase,caseFolder);
%     outputPath = fullfile(outputBase,caseFolder);

%     if ~exist(outputPath,'dir')
%         mkdir(outputPath)
%     end

%     %% =====================================
%     % DIET LOOP
%     %% =====================================

%     for d = 1:length(dietNames)

%         diet = dietNames{d};

%         fprintf('\nProcessing diet: %s\n',diet)

%         %% LOAD FVA FILE

%         filePath = fullfile(inputPath,['FVA_' diet '.mat']);

%         if ~isfile(filePath)
%             fprintf('Skipping %s (file not found)\n',diet);
%             continue;
%         end

%         data = load(filePath);
%         rxns = data.rxns;

%         min_control = data.min_control;
%         max_control = data.max_control;

%         min_PC = data.min_PC;
%         max_PC = data.max_PC;

%         min_PC_opt = data.min_PC_opt;
%         max_PC_opt = data.max_PC_opt;

%         %% =====================================
%         % FIX: MATCH LENGTHS
%         %% =====================================

%         if length(min_control) ~= length(min_PC) || ...
%                 length(min_PC) ~= length(min_PC_opt) || ...
%                 length(rxns) ~= length(min_control)

%             error('Mismatch in FVA output lengths');
%         end

%         n_common = min([length(min_control), length(min_PC), length(min_PC_opt)]);

%         min_control = min_control(1:n_common);
%         max_control = max_control(1:n_common);

%         min_PC = min_PC(1:n_common);
%         max_PC = max_PC(1:n_common);

%         min_PC_opt = min_PC_opt(1:n_common);
%         max_PC_opt = max_PC_opt(1:n_common);

%         if length(rxns) < n_common
%             error('Reaction list shorter than flux vectors');
%         end

%         Reaction = rxns(1:n_common);

%         %% CREATE TABLES

%         fva_CT = table(Reaction, min_control, max_control, ...
%             'VariableNames', {'Reaction','MinFlux','MaxFlux'});

%         fva_PC_tbl = table(Reaction, min_PC, max_PC, ...
%             'VariableNames', {'Reaction','MinFlux','MaxFlux'});

%         fva_PC_opt_tbl = table(Reaction, min_PC_opt, max_PC_opt, ...
%             'VariableNames', {'Reaction','MinFlux','MaxFlux'});

%         %% COMPARISONS

%         comparisons = {
%             fva_PC_tbl,     ['CT_vs_PC_' diet];
%             fva_PC_opt_tbl, ['CT_vs_PCOPT_' diet];
%         };

%         %% LOOP COMPARISONS

%         for c = 1:size(comparisons,1)

%             fva_PC_current = comparisons{c,1};
%             outputTag = comparisons{c,2};

%             fprintf('Running FSR: %s\n',outputTag)

%             commonrxns = fva_CT.Reaction;

%             maxFlux_normal = fva_CT.MaxFlux;
%             minFlux_normal = fva_CT.MinFlux;

%             maxFlux_disease = fva_PC_current.MaxFlux;
%             minFlux_disease = fva_PC_current.MinFlux;

%             %% FSR CALCULATION

%             fluxspanratio = zeros(length(commonrxns),1);

%             for i = 1:length(commonrxns)

%                 if maxFlux_disease(i) ~= maxFlux_normal(i) && ...
%                         minFlux_disease(i) <= minFlux_normal(i)

%                     denom = (maxFlux_normal(i) - minFlux_normal(i));

%                     if denom ~= 0
%                         fluxspanratio(i) = ...
%                             (maxFlux_disease(i) - minFlux_disease(i)) / denom;
%                     else
%                         fluxspanratio(i) = 0;
%                     end

%                 elseif maxFlux_disease(i) == maxFlux_normal(i) && ...
%                         minFlux_disease(i) == minFlux_normal(i)

%                     fluxspanratio(i) = 0;
%                 end
%             end

%             %% CLASSIFICATION

%             Status = repmat("Unchanged", length(commonrxns), 1);

%             up_mask = (fluxspanratio ~= 0 & isfinite(fluxspanratio) & ...
%                 ~isnan(fluxspanratio) & fluxspanratio >= 2 & ...
%                 fluxspanratio < 20000);

%             down_mask = (fluxspanratio ~= 0 & isfinite(fluxspanratio) & ...
%                 ~isnan(fluxspanratio) & fluxspanratio <= 0.5 & ...
%                 fluxspanratio >= 0.01);

%             Status(up_mask) = "Upregulated";
%             Status(down_mask) = "Downregulated";

%             %% STATS

%             num_up = sum(up_mask);
%             num_down = sum(down_mask);
%             total_rxns = length(commonrxns);

%             percent_up = (num_up / total_rxns) * 100;
%             percent_down = (num_down / total_rxns) * 100;

%             fprintf('Total reactions: %d\n', total_rxns);
%             fprintf('Upregulated: %d (%.2f%%)\n', num_up, percent_up);
%             fprintf('Downregulated: %d (%.2f%%)\n', num_down, percent_down);

%             %% TABLE

%             FSR_table = table(commonrxns, fluxspanratio, Status, ...
%                 'VariableNames', {'Reaction','FSR','RegulationStatus'});

%             %% SAVE FILES

%             writetable(FSR_table, ...
%                 fullfile(outputPath,['FSR_' outputTag '.csv']));

%             writetable(FSR_table(Status~="Unchanged",:), ...
%                 fullfile(outputPath,['FSR_Dysregulated_' outputTag '.csv']));

%             writetable(FSR_table, ...
%                 fullfile(outputPath,['FSR_All_' outputTag '.csv']));

%             fprintf('Saved → %s\n',outputTag);

%         end
%     end
% end

% fprintf('\nALL FSR RESULTS SAVED IN "results" FOLDER.\n');

%% =========================================
% DIET LIST
%% =========================================

dietNames = {
    'highfiber'
};

%% =========================================
% INPUT + OUTPUT FOLDERS
%% =========================================

inputBase  = 'FVA_results-HF-allCases';
outputBase = 'FSR_results-HF-allCases';

if ~exist(outputBase,'dir')
    mkdir(outputBase)
end

%% =========================================
% CASE NAMES
%% =========================================

caseNames = {
    'case1_equal'
    'case2_50_50'
    'case3_60_40'
    'case4_70_30'
    'case5_pcos'
};

%% =========================================
% LOOP OVER CASES
%% =========================================

for cs = 1:length(caseNames)

    caseName = caseNames{cs};

    fprintf('\n===========================\n');
    fprintf('Processing %s\n', caseName);
    fprintf('===========================\n');

    inputPath  = fullfile(inputBase, caseName);
    outputPath = fullfile(outputBase, caseName);

    if ~exist(outputPath,'dir')
        mkdir(outputPath)
    end

    %% =====================================
    % DIET LOOP
    %% =====================================

    for d = 1:length(dietNames)

        diet = dietNames{d};
        fprintf('\nProcessing diet: %s\n', diet)

        %% LOAD FVA FILE

        filePath = fullfile(inputPath, ['FVA_' diet '.mat']);

        if ~isfile(filePath)
            fprintf('Skipping %s (%s not found)\n', caseName, diet);
            continue;
        end

        data = load(filePath);

        %% =========================
        % CT vs PC
        %% =========================

        commonrxns = data.commonRxns;
        rxns = string(commonrxns(:));

        min_CT = data.min_CT_common;
        max_CT = data.max_CT_common;

        min_PC = data.min_PC_common;
        max_PC = data.max_PC_common;

        %% =========================
        % CT vs PC_opt
        %% =========================

        commonrxns_opt = data.commonRxns_opt;
        rxns_opt = string(commonrxns_opt(:));

        min_CT_opt = data.min_CT_common_opt;
        max_CT_opt = data.max_CT_common_opt;

        min_PCopt = data.min_PCopt_common;
        max_PCopt = data.max_PCopt_common;

        %% =====================================
        % RUN FSR FOR BOTH COMPARISONS
        %% =====================================

        comparisons = {
            rxns,     min_CT, max_CT, min_PC, max_PC, ['CT_vs_PC_' diet];
            rxns_opt, min_CT_opt, max_CT_opt, min_PCopt, max_PCopt, ['CT_vs_PCOPT_' diet];
        };

        for c = 1:size(comparisons,1)

            rxnList = comparisons{c,1};
            min_normal = comparisons{c,2};
            max_normal = comparisons{c,3};
            min_disease = comparisons{c,4};
            max_disease = comparisons{c,5};
            tag = comparisons{c,6};

            fprintf('Running FSR: %s\n', tag)

            %% =========================
            % SPAN CALCULATIONS
            %% =========================

            nRxns = length(rxnList);

            CT_span = max_normal - min_normal;
            Disease_span = max_disease - min_disease;
            SpanDiff = Disease_span - CT_span;

            fluxspanratio = zeros(nRxns,1);
            AbsDistToCT = zeros(nRxns,1);

            tol = 1e-9;

            for i = 1:nRxns
                denom = CT_span(i);

                if abs(denom) > tol
                    fluxspanratio(i) = Disease_span(i) / denom;
                    AbsDistToCT(i) = abs(fluxspanratio(i) - 1);
                else
                    fluxspanratio(i) = NaN;
                    AbsDistToCT(i) = NaN;
                end
            end

            %% =========================
            % DISTANCE SUMMARY
            %% =========================



            %% =========================
            % TOP 10 CHANGES
            %% =========================

            [sortedDist, sortIdx] = sort(AbsDistToCT, 'descend', 'MissingPlacement', 'last');

            fprintf('Top 10 farthest reactions for %s:\n', tag);

            count = 0;
            for k = 1:length(sortIdx)
                idx = sortIdx(k);

                if ~isnan(AbsDistToCT(idx))
                    fprintf('  %s : %.4f\n', rxnList(idx), AbsDistToCT(idx));
                    count = count + 1;
                end

                if count == 10
                    break;
                end
            end

            %% =========================
            % CLASSIFICATION
            %% =========================

            Status = repmat("Unchanged", nRxns, 1);

            up_mask   = (fluxspanratio >= 2);
            down_mask = (fluxspanratio <= 0.5);

            Status(up_mask)   = "Upregulated";
            Status(down_mask) = "Downregulated";

            Direction = repmat("Unchanged", nRxns, 1);
            Direction(fluxspanratio > 1) = "Up";
            Direction(fluxspanratio < 1) = "Down";

            %% =========================
            % STATS
            %% =========================

            num_up = sum(up_mask, 'omitnan');
            num_down = sum(down_mask, 'omitnan');

            fprintf('Total: %d | Up: %d | Down: %d\n', ...
                nRxns, num_up, num_down);

            %% =========================
            % TABLE
            %% =========================

            FSR_table = table(rxnList, CT_span, Disease_span, SpanDiff, ...
                fluxspanratio, AbsDistToCT, Direction, Status, ...
                'VariableNames', {'Reaction','CT_Span','Disease_Span','SpanDiff', ...
                'FSR','AbsDistToCT','Direction','RegulationStatus'});

            %% =========================
            % SAVE
            %% =========================

            writetable(FSR_table, ...
                fullfile(outputPath, ['FSR_' tag '.csv']));

            writetable(FSR_table(Status~="Unchanged",:), ...
                fullfile(outputPath, ['FSR_Dysregulated_' tag '.csv']));
        end
    end
end

fprintf('\nALL FSR RESULTS SAVED.\n');