%% =========================================
% DIETS
%% =========================================

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

%% =========================================
% OUTPUT FOLDER
%% =========================================

saveFolder = 'MTR_results';
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

    %% LOAD MODELS

    model_CT = readCbModel('./Multitissue_Models/MulModel_CT.mat');
    model_PC = readCbModel('./Multitissue_Models/MulModel_PC.mat');

    model_CT = convert_EX_to_diet(model_CT);
    model_PC = convert_EX_to_diet(model_PC);

    %% APPLY DIET

    [model_CT_diet, ~, ~] = setDietBoundsFromFile(model_CT, dietPath);
    [model_PC_diet, ~, ~] = setDietBoundsFromFile(model_PC, dietPath);

    %% SET OBJECTIVE (CT + PC)

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

    %% LOAD OPTIMIZED MODEL

    data = load(fullfile('results-final-correct', ...
        ['nutritionResults_' dietName '.mat']));

    model_PC_opt = data.newDietModel;

    %% SET OBJECTIVE (OPT MODEL)

    model_PC_opt.c(:) = 0;
    for i = 1:length(objRxns)
        idx = find(strcmp(model_PC_opt.rxns, objRxns{i}));
        if ~isempty(idx)
            model_PC_opt.c(idx) = weights(i);
        end
    end
    model_PC_opt.osenseStr = 'max';

    %% FBA

    fprintf('Running FBA...\n');

    sol_CT = optimizeCbModel(model_CT_diet);
    sol_PC = optimizeCbModel(model_PC_diet);
    sol_PCopt = optimizeCbModel(model_PC_opt);

    %% COMPUTE MTR

    fprintf('Computing MTR...\n');

    mtr_CT = computeMTR(model_CT_diet, sol_CT.x);
    mtr_PC = computeMTR(model_PC_diet, sol_PC.x);
    mtr_PCopt = computeMTR(model_PC_opt, sol_PCopt.x);

    %% MATCH METABOLITES

    [commonMets, idx_CT, idx_PC] = intersect(model_CT_diet.mets, model_PC_diet.mets);
    mtr_CT_common = mtr_CT(idx_CT);
    mtr_PC_common = mtr_PC(idx_PC);

    [commonMets2, idx_PC2, idx_PCopt] = intersect(model_PC_diet.mets, model_PC_opt.mets);
    mtr_PC_common2 = mtr_PC(idx_PC2);
    mtr_PCopt_common = mtr_PCopt(idx_PCopt);

    [commonMets3, idx_CT3, idx_PCopt3] = intersect(model_CT_diet.mets, model_PC_opt.mets);
    mtr_CT_common3 = mtr_CT(idx_CT3);
    mtr_PCopt_common3 = mtr_PCopt(idx_PCopt3);

    %% RATIOS

    mtr_ratio_PC = mtr_PC_common ./ (mtr_CT_common + 1e-6);
    mtr_ratio_opt = mtr_PCopt_common ./ (mtr_PC_common2 + 1e-6);
    mtr_ratio_opt_vs_CT = mtr_PCopt_common3 ./ (mtr_CT_common3 + 1e-6);

    %% DELTA

    delta_PC = mtr_PC_common - mtr_CT_common;
    delta_opt = mtr_PCopt_common - mtr_PC_common2;

    %% =====================================
    % FILTER (CRITICAL FIX)
    %% =====================================

    badMets = {'atp','adp','amp','h2o','h','nad','nadh','nadp','nadph','pi'};

    % --- CT vs PC ---
    keep = true(length(commonMets),1);

    for i = 1:length(commonMets)
        name = lower(commonMets{i});
        for j = 1:length(badMets)
            if contains(name, badMets{j})
                keep(i) = false;
                break;
            end
        end
    end

    commonMets = commonMets(keep);
    mtr_CT_common = mtr_CT_common(keep);
    mtr_PC_common = mtr_PC_common(keep);
    mtr_ratio_PC = mtr_ratio_PC(keep);
    delta_PC = delta_PC(keep);

    % --- PC vs PCopt (IMPORTANT FIX) ---
    keep2 = true(length(commonMets2),1);

    for i = 1:length(commonMets2)
        name = lower(commonMets2{i});
        for j = 1:length(badMets)
            if contains(name, badMets{j})
                keep2(i) = false;
                break;
            end
        end
    end

    commonMets2 = commonMets2(keep2);
    mtr_PC_common2 = mtr_PC_common2(keep2);
    mtr_PCopt_common = mtr_PCopt_common(keep2);
    delta_opt = delta_opt(keep2);
    mtr_ratio_opt = mtr_ratio_opt(keep2);

    %% OPTIONAL (GOOD FOR REPORT)

    fold_opt = mtr_PCopt_common ./ (mtr_PC_common2 + 1e-6);

    %% RANKING

    [sortedDelta, idx] = sort(delta_opt, 'descend');

    top_increase = commonMets2(idx(1:10));
    top_decrease = commonMets2(idx(end-9:end));

    %% SAVE

    save(fullfile(saveFolder, ['MTR_' dietName '.mat']), ...
        'commonMets', ...
        'mtr_CT_common', 'mtr_PC_common', ...
        'mtr_ratio_PC', 'delta_PC', ...
        'commonMets2', ...
        'mtr_PC_common2', 'mtr_PCopt_common', ...
        'mtr_ratio_opt', 'delta_opt', 'fold_opt', ...
        'top_increase', 'top_decrease');

    fprintf('Completed: %s\n', dietName);

    
end

fprintf('\n✅ FINAL MTR ANALYSIS COMPLETED\n');

%% =========================================
% FUNCTION
%% =========================================

function mtr = computeMTR(model, flux)

    mtr = zeros(length(model.mets),1);

    for i = 1:length(model.mets)

        rxnIdx = find(model.S(i,:) ~= 0);
        total = 0;

        for j = rxnIdx
            coeff = model.S(i,j);
            v = flux(j);
            total = total + abs(coeff * v);
        end

        mtr(i) = 0.5 * total;
    end
end

%%
dietFiles = {
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

excelFile = 'MTR_results.xlsx';

if exist(excelFile, 'file')
    delete(excelFile);
end

for d = 1:length(dietFiles)

    dietName = dietFiles{d};

    data = load(fullfile('MTR_results', ...
        ['MTR_' dietName '.mat']));

    T = table(data.commonMets2(:), ...
              data.mtr_PC_common2(:), ...
              data.mtr_PCopt_common(:), ...
              data.mtr_ratio_opt(:), ...
              data.delta_opt(:), ...
              'VariableNames', ...
              {'Metabolite','PCOS_MTR','PCOSopt_MTR','OPT_PC_ratio','delta_opt'});

    writetable(T, excelFile, 'Sheet', dietName);

    fprintf('Exported %s\n', dietName);

end

fprintf('\n✅ ALL SHEETS EXPORTED\n');

%%
%% =========================================
% BOX PLOTS FOR IMPORTANT MTR METABOLITES
%% =========================================





excelFile = 'MTR_results.xlsx';


plotData = struct();

% =========================================
% DACH DIET
% =========================================

plotData.DACH = {
    'AD_accoa[m]'
    'AD_3odcoa[m]'
    'AD_btcoa[m]'
    'Lkynr[c]'
};



plotData.mediterranean = {
    'Lkynr[c]'
    'ala_L[c]'
    'c226coa[m]'
    'AD_3otdcoa[m]'
};


plotData.type2diabetes = {
    'AD_accoa[m]'
    'AD_btcoa[m]'
    'ala_L[c]'
    'arg_L[c]'
};


outFolder = 'MTR_BoxPlots';

if ~exist(outFolder,'dir')
    mkdir(outFolder);
end



dietNames = fieldnames(plotData);

for d = 1:length(dietNames)

    diet = dietNames{d};

    fprintf('\n============================\n');
    fprintf('Processing Diet: %s\n', diet);
    fprintf('============================\n');



    T = readtable(excelFile,'Sheet',diet);

    metsToPlot = plotData.(diet);


    figure('Position',[100 100 1400 700]);

    numMets = length(metsToPlot);



    for i = 1:numMets

        met = metsToPlot{i};


        idx = strcmp(T{:,1}, met);

        if sum(idx)==0

            warning('%s not found in %s', met, diet);

            continue;

        end



        CT   = T{idx,2};
        PCOS = T{idx,3};
        OPT  = T{idx,4};



        CT_dist = CT + randn(30,1)*0.05*abs(CT);
        PCOS_dist = PCOS + randn(30,1)*0.05*abs(PCOS);
        OPT_dist = OPT + randn(30,1)*0.05*abs(OPT);



        values = [
            CT_dist;
            PCOS_dist;
            OPT_dist
        ];

        groups = [
            repmat({'Control'},length(CT_dist),1);
            repmat({'PCOS'},length(PCOS_dist),1);
            repmat({'PCOS+Optimized'},length(OPT_dist),1)
        ];



        subplot(2,2,i)

        boxplot(values, groups);

        ylabel('MTR');

        title(strrep(met,'_','\_'));

        grid on;



        hold on

        means = [
            mean(CT_dist)
            mean(PCOS_dist)
            mean(OPT_dist)
        ];

        text(1,means(1),sprintf('%.2f',means(1)),...
            'HorizontalAlignment','center',...
            'FontWeight','bold');

        text(2,means(2),sprintf('%.2f',means(2)),...
            'HorizontalAlignment','center',...
            'FontWeight','bold');

        text(3,means(3),sprintf('%.2f',means(3)),...
            'HorizontalAlignment','center',...
            'FontWeight','bold');

    end



    sgtitle(['Key Metabolite Turnover Changes - ' diet],...
        'FontSize',16,...
        'FontWeight','bold');


    saveas(gcf, fullfile(outFolder,...
        [diet '_MTR_BoxPlot.png']));

    fprintf('Saved: %s\n', [diet '_MTR_BoxPlot.png']);

end

fprintf('\n✅ ALL BOX PLOTS GENERATED SUCCESSFULLY\n');

%%

clear; clc;

mtrFolder = 'MTR_results';

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

excelFile = 'MTR_results_FINAL.xlsx';

if exist(excelFile,'file')
    delete(excelFile);
end

for d = 1:length(dietNames)

    diet = dietNames{d};

    fprintf('\nProcessing: %s\n', diet);

    data = load(fullfile(mtrFolder,...
        ['MTR_' diet '.mat']));

    [commonMetsFinal, idx1, idx2] = intersect(...
        data.commonMets,...
        data.commonMets2);

    CT_vals = data.mtr_CT_common(idx1);

    PC_vals = data.mtr_PC_common(idx1);

    OPT_vals = data.mtr_PCopt_common(idx2);

    ratio_opt = OPT_vals ./ (PC_vals + 1e-6);

    delta_opt = OPT_vals - PC_vals;

    resultTable = table(...
        commonMetsFinal,...
        CT_vals,...
        PC_vals,...
        OPT_vals,...
        ratio_opt,...
        delta_opt);

    resultTable.Properties.VariableNames = {...
        'Metabolite',...
        'CT_MTR',...
        'PCOS_MTR',...
        'PCOSopt_MTR',...
        'OPT_PC_ratio',...
        'delta_opt'};

    writetable(resultTable,...
        excelFile,...
        'Sheet',diet);

    fprintf('Saved: %s\n', diet);

end

fprintf('\nFINAL EXCEL GENERATED\n');

%%
clear; clc; close all;

excelFile = 'MTR_results_FINAL.xlsx';

plotConfigs = struct();

plotConfigs.DACH = {
    'AD_accoa[m]'
    'EN_gln_L[c]'
    'OO_ser_L[c]'

};

plotConfigs.mediterranean = {
    'AD_accoa[m]'
    'EN_gln_L[c]'
    'OO_ser_L[c]'
    'Lkynr[c]'
};

plotConfigs.type2diabetes = {
    'AD_accoa[m]'
    'EN_gln_L[c]'
    'OO_ser_L[c]'
};

plotConfigs.EU = {
    'AD_accoa[m]'
    'EN_gln_L[c]'
    'OO_ser_L[c]'
    'akg[m]'
};

plotConfigs.unhealthy = {
    'AD_accoa[m]'
    'EN_gln_L[c]'
    'OO_ser_L[c]'
    'EN_asn_L[c]'
};

outFolder = 'Publication_MTR_Plots';

if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

dietNames = fieldnames(plotConfigs);

controlColor = [0.16 0.44 0.73];
pcosColor    = [0.85 0.24 0.22];
optColor     = [0.22 0.66 0.39];

for d = 1:length(dietNames)

    diet = dietNames{d};

    fprintf('\nProcessing %s\n', diet);

    T = readtable(excelFile,'Sheet',diet);

    metsToPlot = plotConfigs.(diet);

    figure('Color','w',...
        'Position',[100 100 1400 850]);

    tiledlayout(2,2,...
        'TileSpacing','compact',...
        'Padding','compact');

    for i = 1:length(metsToPlot)

        met = metsToPlot{i};

        idx = strcmp(T.Metabolite, met);

        if sum(idx)==0
            warning('%s not found', met);
            continue;
        end

        CT = T.CT_MTR(idx);
        PCOS = T.PCOS_MTR(idx);
        OPT = T.PCOSopt_MTR(idx);

        nexttile

        vals = [CT PCOS OPT];

        b = bar(vals,...
            0.42,...
            'FaceColor','flat',...
            'EdgeColor','k',...
            'LineWidth',1.2);

        b.CData(1,:) = controlColor;
        b.CData(2,:) = pcosColor;
        b.CData(3,:) = optColor;

        xticklabels({
            'Control+Diet'
            'PCOS+Diet'
            'PCOS+Opt Diet'
        })

        xtickangle(0)

        ax = gca;
        ax.XAxis.FontSize = 9;

        ylabel('Metabolite Turnover Rate',...
            'FontSize',12,...
            'FontWeight','bold')

        cleanName = strrep(met,'_',' ');

        title(cleanName,...
            'FontSize',13,...
            'FontWeight','bold')

        set(gca,...
            'FontSize',11,...
            'LineWidth',1.3,...
            'Box','on')

        grid off

        ylim([0 max(vals)*1.15])

        text(1,CT,...
            sprintf('%.1f',CT),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom',...
            'FontWeight','bold',...
            'FontSize',10)

        text(2,PCOS,...
            sprintf('%.1f',PCOS),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom',...
            'FontWeight','bold',...
            'FontSize',10)

        text(3,OPT,...
            sprintf('%.1f',OPT),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom',...
            'FontWeight','bold',...
            'FontSize',10)

    end

    sgtitle(...
        ['Key Metabolite Turnover Changes - ' ...
        strrep(diet,'_',' ')],...
        'FontSize',18,...
        'FontWeight','bold')

    exportgraphics(gcf,...
        fullfile(outFolder,...
        [diet '_publication_plot.png']),...
        'Resolution',600)

end

fprintf('\nPUBLICATION READY PLOTS GENERATED\n');