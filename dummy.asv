%% ===========================
% 1. LOAD MODEL
% ===========================
modelData = load('Individual_models/AD_CT.mat');
fn = fieldnames(modelData);
model = modelData.(fn{1});

%% ===========================
% 2. EXTRACT SUBSYSTEMS
% ===========================
rawSubs = model.subSystems;
subs = cell(size(rawSubs));
for i = 1:length(rawSubs)
    if iscell(rawSubs{i})
        subs{i} = rawSubs{i}{1};
    else
        subs{i} = rawSubs{i};
    end
end
subs = subs(~cellfun('isempty', subs));

%% ===========================
% 3. SELECT FATTY-ACID SUBSYSTEMS
% ===========================
targetSubsystems = {"Fatty acid oxidation", "Fatty acid synthesis"};
selectedRxns = {};

for t = 1:length(targetSubsystems)
    ss = targetSubsystems{t};
    idx = find(strcmp(subs, ss));
    selectedRxns = [selectedRxns; model.rxns(idx)];
end

selectedRxns = unique(selectedRxns);
fprintf('Total selected fatty-acid reactions: %d\n', length(selectedRxns));

%% ===========================
% 4. DEFINE FATTY-ACID OBJECTIVE
% ===========================
modelFA = model;
modelFA.c(:) = 0;

for i = 1:length(selectedRxns)
    idx_r = find(strcmp(modelFA.rxns, selectedRxns{i}));
    if ~isempty(idx_r)
        modelFA.c(idx_r) = 1;
    end
end

modelFA.osenseStr = 'max';
solFA = optimizeCbModel(modelFA);
fprintf('\nFatty acid objective flux (WT) = %.6f\n', solFA.f);

%% ===========================
% 5. FVA UNDER NEW OBJECTIVE
% ===========================
[minV, maxV] = fluxVariability(modelFA, 100, 'max', selectedRxns);
active = find(maxV > 1e-6);
selectedRxns = selectedRxns(active);

fprintf('Active reactions after FVA filtering: %d\n', length(selectedRxns));

%% ===========================
% 6. SINGLE REACTION DELETION
% ===========================
fprintf('\nRunning single reaction deletion...\n');

[out1, out2] = singleRxnDeletion(modelFA, 'FBA', selectedRxns);

if isstruct(out1)
    grWT = out1.grRateWT;
    grKO = out1.grRateKO;
else
    grKO = out1;
    grWT = out2;
end

essScore = grKO ./ grWT;
importantIdx = find(essScore < 0.5);
importantReactions = selectedRxns(importantIdx);

fprintf('Reactions significantly affecting FA objective (>50%% drop): %d\n', length(importantReactions));

%% Stop if nothing important
if isempty(importantReactions)
    fprintf('\nNo important reactions found. Exiting.\n');
    return
end

%% ===========================
% 7. RANK IMPORTANT REACTIONS
% ===========================
[minImp, maxImp] = fluxVariability(modelFA, 100, 'max', importantReactions);
fluxSpan = maxImp - minImp;

ess_imp = essScore(importantIdx);
flux_imp = fluxSpan(:);

n = min([length(importantReactions), length(ess_imp), length(flux_imp)]);
importantReactions = importantReactions(1:n);
ess_imp = ess_imp(1:n);
flux_imp = flux_imp(1:n);

importanceScore = (1 - ess_imp) .* (1 ./ (flux_imp + 1e-9));

[~, sidx] = sort(importanceScore, 'descend');
topReactions = importantReactions(sidx(1));

fprintf('\n===== TOP CONTROLLING REACTIONS IN FAT METABOLISM =====\n');
disp(topReactions);

%% ===========================
% 8. SAVE FULL RESULTS
% ===========================
T = table(importantReactions, ess_imp, flux_imp, importanceScore, ...
    'VariableNames', {'Reaction', 'EssentialityDrop', 'FluxSpan', 'ImportanceScore'});

writetable(T, 'FattyAcid_KeyReactions.csv');
fprintf('\nSaved detailed results to FattyAcid_KeyReactions.csv\n');
