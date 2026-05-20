%% =========================================
% DIETS
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

inputFolder  = 'FVA_results_ROI';
outputFolder = 'FSR_results_ROI';

if ~exist(outputFolder,'dir')
    mkdir(outputFolder)
end

%% =========================================
% LOAD BASELINE
%% =========================================

base = load(fullfile(inputFolder,'FVA_ROI_BASELINE.mat'));

roi_base = string(base.roi_valid(:));

span_CT_base = base.span_CT_base;
span_PC_base = base.span_PC_base;

%% =========================================
% LOOP OVER DIETS
%% =========================================

for d = 1:length(dietNames)

    diet = dietNames{d};
    fprintf('\n===========================\n');
    fprintf('Processing diet: %s\n', diet);
    fprintf('===========================\n');

    filePath = fullfile(inputFolder, ['FVA_ROI_' diet '.mat']);

    if ~isfile(filePath)
        fprintf('Skipping %s (file not found)\n', diet);
        continue;
    end

    data = load(filePath);

    rxns = string(data.roi_valid(:));

    span_CT    = data.span_CT;
    span_PC    = data.span_PC;
    span_PCopt = data.span_PCopt;

    %% =====================================
    % ALIGN WITH BASELINE
    %% =====================================

    [commonRxns, idx_base, idx_curr] = intersect(roi_base, rxns, 'stable');

    span_CT_base_m = span_CT_base(idx_base);
    span_PC_base_m = span_PC_base(idx_base);

    span_CT_m    = span_CT(idx_curr);
    span_PC_m    = span_PC(idx_curr);
    span_PCopt_m = span_PCopt(idx_curr);

    %% =====================================
    % FSR CALCULATION
    %% =====================================

    tol = 1e-9;

    FSR_PC = span_PC_m ./ span_CT_m;
    FSR_PC(abs(span_CT_m) <= tol) = NaN;

    FSR_PCOPT = span_PCopt_m ./ span_CT_m;
    FSR_PCOPT(abs(span_CT_m) <= tol) = NaN;

    FSR_BASE_PC = span_PC_m ./ span_CT_base_m;
    FSR_BASE_PC(abs(span_CT_base_m) <= tol) = NaN;

    FSR_BASE_PCOPT = span_PCopt_m ./ span_CT_base_m;
    FSR_BASE_PCOPT(abs(span_CT_base_m) <= tol) = NaN;

    %% =====================================
    % CLASSIFICATION (NO FUNCTION)
    %% =====================================

    % CT vs PC
    Status_PC = repmat("Unchanged", length(FSR_PC), 1);
    Status_PC(FSR_PC >= 2) = "Upregulated";
    Status_PC(FSR_PC <= 0.5 & FSR_PC >= 0.01) = "Downregulated";

    % CT vs PCOPT
    Status_PCOPT = repmat("Unchanged", length(FSR_PCOPT), 1);
    Status_PCOPT(FSR_PCOPT >= 2) = "Upregulated";
    Status_PCOPT(FSR_PCOPT <= 0.5 & FSR_PCOPT >= 0.01) = "Downregulated";

    % BASE vs PC
    Status_BASE_PC = repmat("Unchanged", length(FSR_BASE_PC), 1);
    Status_BASE_PC(FSR_BASE_PC >= 2) = "Upregulated";
    Status_BASE_PC(FSR_BASE_PC <= 0.5 & FSR_BASE_PC >= 0.01) = "Downregulated";

    % BASE vs PCOPT
    Status_BASE_PCOPT = repmat("Unchanged", length(FSR_BASE_PCOPT), 1);
    Status_BASE_PCOPT(FSR_BASE_PCOPT >= 2) = "Upregulated";
    Status_BASE_PCOPT(FSR_BASE_PCOPT <= 0.5 & FSR_BASE_PCOPT >= 0.01) = "Downregulated";

    %% =====================================
    % STATS + PERCENTAGES
    %% =====================================

    n = length(commonRxns);

    % ---- CT vs PC ----
    up = sum(FSR_PC >= 2, 'omitnan');
    down = sum(FSR_PC <= 0.5 & FSR_PC >= 0.01, 'omitnan');

    fprintf('CT vs PC        → Up:%d (%.2f%%) | Down:%d (%.2f%%) | Total: %.2f%%\n', ...
        up, (up/n)*100, down, (down/n)*100, ((up+down)/n)*100);

    % ---- CT vs PCOPT ----
    up = sum(FSR_PCOPT >= 2, 'omitnan');
    down = sum(FSR_PCOPT <= 0.5 & FSR_PCOPT >= 0.01, 'omitnan');

    fprintf('CT vs PCOPT     → Up:%d (%.2f%%) | Down:%d (%.2f%%) | Total: %.2f%%\n', ...
        up, (up/n)*100, down, (down/n)*100, ((up+down)/n)*100);

    % ---- BASE vs PC ----
    up = sum(FSR_BASE_PC >= 2, 'omitnan');
    down = sum(FSR_BASE_PC <= 0.5 & FSR_BASE_PC >= 0.01, 'omitnan');

    fprintf('BASE vs PC      → Up:%d (%.2f%%) | Down:%d (%.2f%%) | Total: %.2f%%\n', ...
        up, (up/n)*100, down, (down/n)*100, ((up+down)/n)*100);

    % ---- BASE vs PCOPT ----
    up = sum(FSR_BASE_PCOPT >= 2, 'omitnan');
    down = sum(FSR_BASE_PCOPT <= 0.5 & FSR_BASE_PCOPT >= 0.01, 'omitnan');

    fprintf('BASE vs PCOPT   → Up:%d (%.2f%%) | Down:%d (%.2f%%) | Total: %.2f%%\n', ...
        up, (up/n)*100, down, (down/n)*100, ((up+down)/n)*100);

    %% =====================================
    % SAVE TABLES
    %% =====================================

    T1 = table(commonRxns, FSR_PC, Status_PC);
    T2 = table(commonRxns, FSR_PCOPT, Status_PCOPT);
    T3 = table(commonRxns, FSR_BASE_PC, Status_BASE_PC);
    T4 = table(commonRxns, FSR_BASE_PCOPT, Status_BASE_PCOPT);

    writetable(T1, fullfile(outputFolder,['FSR_CT_vs_PC_' diet '.csv']));
    writetable(T2, fullfile(outputFolder,['FSR_CT_vs_PCOPT_' diet '.csv']));
    writetable(T3, fullfile(outputFolder,['FSR_BASE_vs_PC_' diet '.csv']));
    writetable(T4, fullfile(outputFolder,['FSR_BASE_vs_PCOPT_' diet '.csv']));

end

fprintf('\n✅ FSR COMPLETED SUCCESSFULLY\n');

%%
clear; clc;

dietNames = {'High Protein','Type2 Diabetes','Mediterranean','DACH','High Fiber',...
    'Vegan','Unhealthy','High Fat Low Carb','Gluten Free','Vegetarian','EU'};

before = [42.40,43.24,42.90,42.57,42.40,42.24,44.74,44.41,43.41,41.24,43.41];
after  = [36.39,37.80,37.73,37.70,37.85,37.90,40.57,40.90,40.40,39.40,42.07];

figure('Color','w','Position',[100 100 1100 450])

barData = [before; after]';
b = bar(barData,'grouped','BarWidth',0.75);

% Clean colors (journal style)
b(1).FaceColor = [0.65 0.65 0.65];   % grey
b(2).FaceColor = [0.2 0.45 0.75];    % blue

% X-axis formatting
xticks(1:length(dietNames))
xticklabels(dietNames)
xtickangle(35)

% Labels
ylabel('Dysregulation (%)','FontSize',12)

% Move legend OUTSIDE (no overlap)
legend({'Initial','After'},...
    'Location','northoutside',...
    'Orientation','horizontal',...
    'Box','off')

% Axis styling
set(gca,'FontSize',11,'LineWidth',1.2)
box off
ylim([34 46])   % tighten range → better visual

title('Effect of Dietary Optimization on Metabolic Dysregulation',...
    'FontWeight','bold','FontSize',13)

% Add spacing so nothing overlaps
ax = gca;
ax.Position = [0.06 0.18 0.9 0.65];

exportgraphics(gcf,'Figure1_Clean.png','Resolution',300)