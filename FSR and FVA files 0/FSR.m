

%% Load Control FVA
fva_CT = readtable('FVA_CT_EU.csv');

%% Define cases
cases = {
    'PC_HF',           'FVA_PC_HF.csv',           'CT_vs_PC_HF';
    'PC_OptimizedHF',  'FVA_PC_OptimizedHF.csv',  'CT_vs_PC_OptimizedHF'
};

%% Loop through both cases
for caseIdx = 1:size(cases,1)

    comparisonMode = cases{caseIdx,1};
    fva_PC_file = cases{caseIdx,2};
    outputTag = cases{caseIdx,3};

    fprintf('\n====================================\n');
    fprintf('Running FSR Case: %s\n', comparisonMode);
    fprintf('====================================\n');

    %% Load diseased model FVA
    fva_PC = readtable(fva_PC_file);

    %% Find common reactions
    [common_rxns, idx_CT, idx_PC] = intersect(fva_CT.Reaction, fva_PC.Reaction, 'stable');

    fva_CT_common = fva_CT(idx_CT, :);
    fva_PC_common = fva_PC(idx_PC, :);

    commonrxns = common_rxns;

    maxFlux_normal_met   = fva_CT_common.MaxFlux;
    minFlux_normal_met   = fva_CT_common.MinFlux;

    maxFlux_diseased_met = fva_PC_common.MaxFlux;
    minFlux_diseased_met = fva_PC_common.MinFlux;

    %% ================================
    % ORIGINAL FSR LOGIC — UNCHANGED
    %% ================================

    fluxspanrationew0 = [];

    for i = 1:length(commonrxns)
        if maxFlux_diseased_met(i) ~= maxFlux_normal_met(i) && ...
           minFlux_diseased_met(i) <= minFlux_normal_met(i)

            fluxspanrationew0(i) = ...
                (maxFlux_diseased_met(i) - minFlux_diseased_met(i)) ./ ...
                (maxFlux_normal_met(i) - minFlux_normal_met(i));

        elseif maxFlux_diseased_met(i) == maxFlux_normal_met(i) && ...
               minFlux_diseased_met(i) == minFlux_normal_met(i)

            fluxspanrationew0(i) = 0;
        end
    end

    fluxspanrationew = fluxspanrationew0';

    %% ================================
    % CLASSIFY REGULATION STATUS
    %% ================================

    Downregulated_disease = commonrxns( ...
        find((fluxspanrationew ~= 0 & isfinite(fluxspanrationew) & ...
              ~isnan(fluxspanrationew) & fluxspanrationew <= 0.5 & ...
              fluxspanrationew >= 0.01)));

    Upregulated_disease = commonrxns( ...
        fluxspanrationew ~= 0 & isfinite(fluxspanrationew) & ...
        ~isnan(fluxspanrationew) & fluxspanrationew >= 2 & ...
        fluxspanrationew < 20000.0);

    fsr_Downregulated_disease = fluxspanrationew( ...
        find((fluxspanrationew ~= 0 & isfinite(fluxspanrationew) & ...
              ~isnan(fluxspanrationew) & fluxspanrationew <= 0.5 & ...
              fluxspanrationew >= 0.01)));

    fsr_Upregulated_disease = fluxspanrationew( ...
        find((fluxspanrationew ~= 0 & isfinite(fluxspanrationew) & ...
              ~isnan(fluxspanrationew) & fluxspanrationew >= 2 & ...
              fluxspanrationew < 20000.0)));

    %% ================================
    % COUNT & PERCENT
    %% ================================

    num_up   = length(Upregulated_disease);
    num_down = length(Downregulated_disease);
    total_rxns = length(commonrxns);

    percent_up   = (num_up / total_rxns) * 100;
    percent_down = (num_down / total_rxns) * 100;

    fprintf('Total reactions: %d\n', total_rxns);
    fprintf('Upregulated: %d (%.2f%%)\n', num_up, percent_up);
    fprintf('Downregulated: %d (%.2f%%)\n', num_down, percent_down);

    %% ================================
    % RESULT TABLE
    %% ================================

    Status = repmat("Unchanged", total_rxns, 1);

    up_mask = (fluxspanrationew ~= 0 & isfinite(fluxspanrationew) & ...
               ~isnan(fluxspanrationew) & fluxspanrationew >= 2 & ...
               fluxspanrationew < 20000.0);

    down_mask = (fluxspanrationew ~= 0 & isfinite(fluxspanrationew) & ...
                 ~isnan(fluxspanrationew) & fluxspanrationew <= 0.5 & ...
                 fluxspanrationew >= 0.01);

    Status(up_mask) = "Upregulated";
    Status(down_mask) = "Downregulated";

    FSR_table = table(commonrxns, fluxspanrationew, Status, ...
        'VariableNames', {'Reaction', 'FSR', 'RegulationStatus'});

    %% ================================
    % SAVE OUTPUTS
    %% ================================

    writetable(FSR_table, ['FSR_results_' outputTag '.csv']);

    dysregulated_mask = (Status ~= "Unchanged");
    FSR_dysregulated = FSR_table(dysregulated_mask, :);

    writetable(FSR_dysregulated, ['FSR_Dysregulated_' outputTag '.csv']);

    FSR_all = FSR_table;
    writetable(FSR_all, ['FSR_All_' outputTag '.csv']);

    fprintf('Saved results for %s\n', outputTag);

end

fprintf('\nALL FSR CASES COMPLETED SUCCESSFULLY.\n');
