
%initCobraToolbox(false);
%changeCobraSolver('gurobi','LP');
 
model_CT = readCbModel('./Multitissue_Models/MulModel_CT.mat');
model_PC = readCbModel('./Multitissue_Models/MulModel_PC.mat');

model_CT = convert_EX_to_diet(model_CT);
model_PC = convert_EX_to_diet(model_PC);


eu_dietPath = fullfile('Diets','EU.tsv');
hf_dietPath = fullfile('Diets','highfiber.tsv')
[model_CT_EU, ~, ~] = setDietBoundsFromFile(model_CT, eu_dietPath);
[model_PC_EU, ~, ~] = setDietBoundsFromFile(model_PC, hf_dietPath);
%%

exc_rxns_idx = find(startsWith(model_CT_EU.rxns, 'Diet_EX_'));
zeroUptakeIdx = exc_rxns_idx(model_CT_EU.lb(exc_rxns_idx) == 0);
zeroUptakeRxns = model_CT_EU.rxns(zeroUptakeIdx);

% disp('Exchange reactions with zero uptake (not allowed by diet):');
% disp(zeroUptakeRxns);

% Weighted tissue biomass maintenance objective
% tissueBiomassRxns = {'SK_biomass_maintenance','AD_biomass_maintenance','GN_biomass_maintenance','OO_biomass_maintenance','EN_biomass_maintenance'};
% weights = [0.2, 0.2, 0.2, 0.2, 0.2];
tissueBiomassRxns = {'SK_ATPtm','AD_ACCOAC','AD_biomass_maintenance','GN_biomass_maintenance','GN_P450SCC1m','OO_biomass_maintenance','EN_biomass_maintenance','EN_biomass_reaction'};
weights = [0.2, 0.1,0.1, 0.1, 0.1, 0.2, 0.1, 0.1]; 
model_CT_EU.c(:) = 0;
for i = 1:length(tissueBiomassRxns)
    rxnIdx = find(strcmp(model_CT_EU.rxns, tissueBiomassRxns{i}));
    if ~isempty(rxnIdx)  
        model_CT_EU.c(rxnIdx) = weights(i);
    else
        fprintf('Warning: %s not found in model_CT_EU.rxns\n', tissueBiomassRxns{i});
    end
end
sol_CT_weighted = optimizeCbModel(model_CT_EU, 'max');
fprintf('Weighted tissue biomass objective value from main file in CT: %.6f\n', sol_CT_weighted.f);


biomassRxns = {'SK_biomass_maintenance', 'AD_biomass_maintenance', ...
               'GN_biomass_maintenance', 'OO_biomass_maintenance', ...
               'EN_biomass_maintenance',  'SK_PGESr', ...
               'GN_PGSr', 'GN_HMR_2581', 'OO_HMR_0987', ...
               'OO_LTC4CP', 'GN_RE1796R', 'EN_HAS2', ...
               'AD_ARGSL', 'AD_SPHK11'};

% Find their indices in the model
[isPresent, rxnIndices] = ismember(biomassRxns, model_CT_EU.rxns);

% Print fluxes for each reaction
fprintf('\nFlux values for key CT model reactions:\n');
for i = 1:length(biomassRxns)
    if isPresent(i)
        fprintf('%s = %.4f\n', biomassRxns{i}, sol_CT_weighted.x(rxnIndices(i)));
    else
        fprintf('%s not found in model_CT_EU\n', biomassRxns{i});
    end
end

% commonRxns = intersect(model_CT_EU.rxns, model_PC_EU.rxns, 'stable');
% isDiet = startsWith(commonRxns,'Diet_EX_');
% isExit = startsWith(commonRxns,'Exit_EX_');
% roi_list = commonRxns(~(isDiet | isExit));

% roi_list = {
%     'SK_PGESr'
%     'GN_PGSr'
%     'GN_HMR_2581'
%     'OO_HMR_0987'
%     'OO_LTC4CP'
%     'GN_RE1796R'
%     'EN_HAS2'
%     'AD_ARGSL'
%     'AD_SPHK11'
% };

% Load ROI list from Excel file
roiTable = readtable('rois.xlsx');
roi_list = roiTable{:,1};                % assumes ROIs are in the first column



control_flux_CT_EU = zeros(length(roi_list),1);
for i = 1:length(roi_list)
    rxnIdx = find(strcmp(model_CT_EU.rxns, roi_list{i}));
    control_flux_CT_EU(i) = sol_CT_weighted.v(rxnIdx);
end

options.display = 'on';
options.roiWeights = ones(1,length(roi_list));
options.targetFlux = control_flux_CT_EU;

% Restrict menu changes to only Diet_EX reactions present in the EU diet file
allowedDietRxns = model_CT_EU.rxns(model_CT_EU.lb < 0 & startsWith(model_CT_EU.rxns, 'Diet_EX_'));
options.targetedDietRxns = [allowedDietRxns, num2cell(ones(length(allowedDietRxns),1))];

[newDietModel, pointsModel, roiFlux, pointsModelSln, menuChanges, detailedAnalysis] = ...
    nutritionAlgorithm(model_PC_EU, roi_list, options);

disp('Suggested dietary changes:');
disp(menuChanges);

save('nutritionResults_PC_EU.mat', ...
     'newDietModel', 'pointsModel', 'roiFlux', ...
     'pointsModelSln', 'menuChanges', 'detailedAnalysis');

menuTable = struct2table(menuChanges);
writetable(menuTable, 'menuChanges_PC_EU.xlsx');

