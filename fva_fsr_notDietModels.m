
%% =========================================
% LOAD MODELS
%% =========================================

model_CT = readCbModel('./Multitissue_Models/MulModel_CT.mat');
model_PC = readCbModel('./Multitissue_Models/MulModel_PC.mat');

%% =========================================
% OBJECTIVE (TISSUE)
%% =========================================

tissueRxns = {
    'SK_ATPtm'
    'AD_ACCOAC'
    'GN_biomass_maintenance'
    'GN_P450SCC1m'
    'OO_ATPtm'
    'EN_SERPT'
    'EN_SMS'
};

tissueWeights = [0.2 0.2 0.1 0.1 0.2 0.1 0.1];

% Reset objectives
model_CT.c(:) = 0;
model_PC.c(:) = 0;

% Assign weighted objective
for i = 1:length(tissueRxns)
    idx_CT = strcmp(model_CT.rxns, tissueRxns{i});
    idx_PC = strcmp(model_PC.rxns, tissueRxns{i});

    model_CT.c(idx_CT) = tissueWeights(i);
    model_PC.c(idx_PC) = tissueWeights(i);
end

model_CT.osenseStr = 'max';
model_PC.osenseStr = 'max';

%% =========================================
% RUN FVA
%% =========================================

fprintf('Running FVA: Control (no diet)\n')
[min_CT, max_CT] = fluxVariability(model_CT, 100);

fprintf('Running FVA: PCOS (no diet)\n')
[min_PC, max_PC] = fluxVariability(model_PC, 100);


%% MATCH LENGTHS (FIX)


n_common = min([length(min_CT), length(min_PC)]);

min_CT = min_CT(1:n_common);
max_CT = max_CT(1:n_common);

min_PC = min_PC(1:n_common);
max_PC = max_PC(1:n_common);

Reaction = (1:n_common)';

% FSR CALCULATION


fluxspanratio = zeros(n_common,1);

for i = 1:n_common

    if max_PC(i) ~= max_CT(i) && min_PC(i) <= min_CT(i)

        denom = (max_CT(i) - min_CT(i));

        if denom ~= 0
            fluxspanratio(i) = ...
                (max_PC(i) - min_PC(i)) / denom;
        else
            fluxspanratio(i) = 0;
        end

    elseif max_PC(i) == max_CT(i) && min_PC(i) == min_CT(i)
        fluxspanratio(i) = 0;
    end
end


% CLASSIFICATION


Status = repmat("Unchanged", n_common, 1);

up_mask = (fluxspanratio >= 2);
down_mask = (fluxspanratio <= 0.5 & fluxspanratio >= 0.01);

Status(up_mask) = "Upregulated";
Status(down_mask) = "Downregulated";


% DYSREGULATION %


num_up = sum(up_mask);
num_down = sum(down_mask);
total = n_common;

percent_up = (num_up / total) * 100;
percent_down = (num_down / total) * 100;
total_dys = percent_up + percent_down;

fprintf('\n===== BASELINE (NO DIET, TISSUE OBJECTIVE) =====\n');
fprintf('Total reactions: %d\n', total);
fprintf('Upregulated: %.2f%%\n', percent_up);
fprintf('Downregulated: %.2f%%\n', percent_down);
fprintf('TOTAL DYSREGULATION: %.2f%%\n', total_dys);


% SAVE RESULTS


FSR_table = table(Reaction, fluxspanratio, Status, ...
    'VariableNames', {'Reaction','FSR','RegulationStatus'});

if ~exist('results','dir')
    mkdir('results')
end

writetable(FSR_table, 'results/FSR_noDiet_tissue.csv');
writetable(FSR_table(Status~="Unchanged",:), ...
    'results/FSR_noDiet_dysregulated_tissue.csv');

save('results/FVA_noDiet_tissue.mat','min_CT','max_CT','min_PC','max_PC');

fprintf('\nFVA + FSR COMPLETED.\n');