function [model, totalDietRxns, unmatchedCount] = setDietBoundsFromFile(model, dietFilePath)
fprintf('\n========== SETTING DIET BOUNDS ==========\n');
fprintf('Loaded diet file: %s\n', dietFilePath);

% read diet file
dietTable = readtable(dietFilePath, 'FileType', 'text', 'Delimiter', '\t', 'VariableNamingRule', 'preserve');
totalDietRxns = height(dietTable);
reactionNames = dietTable.Reaction;
fluxValues = dietTable.("Flux Value");
fprintf(' → Total diet entries in file: %d\n', totalDietRxns);

% ensure Diet_EX_ are present (convert EX_ if needed)
if ~any(startsWith(model.rxns, 'Diet_EX_'))
    if any(startsWith(model.rxns, 'EX_'))
        fprintf('Converting EX_ reactions to Diet_EX_...\n');
        model = convert_EX_to_diet(model);
    end
end

% find all exchange-like reactions and reset lower bounds
exc_rxns_idx = find(startsWith(model.rxns, 'Diet_EX_'));
ex_rxns_idx  = find(startsWith(model.rxns, 'EX_'));
all_ex_idx = unique([exc_rxns_idx; ex_rxns_idx]);
if ~isempty(all_ex_idx)
    model.lb(all_ex_idx) = 0; % zero uptake for all exchanges by default
end

% prevent export via Diet_EX_ and avoid invalid ub<0 after convert
if ~isempty(exc_rxns_idx)
    model.ub(exc_rxns_idx) = 1000;
end

% prepare diet metabolite names (strip EX_ if present in file)
fprintf('\n→ Cleaning diet metabolite names...\n');
dietMets = reactionNames;
for i = 1:numel(dietMets)
    if startsWith(dietMets{i}, 'EX_')
        dietMets{i} = extractAfter(dietMets{i}, 3);
    end
end

% get metabolite names associated with Diet_EX_ reactions (remove Bl_ prefix)
fprintf('→ Extracting model external metabolites...\n');
externalMets = cell(numel(exc_rxns_idx), 1);
for i = 1:numel(exc_rxns_idx)
    rxnIdx = exc_rxns_idx(i);
    metIdx = find(model.S(:, rxnIdx) ~= 0);
    if ~isempty(metIdx)
        externalMets{i}  = erase(model.mets{metIdx(1)}, 'Bl_');


    else
        externalMets{i} = '';
    end
end

% match diet file mets to model Diet_EX_ and set bounds
fprintf('\n→ Matching diet metabolites with model Diet_EX_ reactions...\n');
matched = false(numel(dietMets),1);
for i = 1:numel(dietMets)
    idx = find(strcmp(externalMets, dietMets{i}), 1);
    if ~isempty(idx)
        rxnIdx = exc_rxns_idx(idx);
        model.lb(rxnIdx) = -(fluxValues(i)/(24*10000));
        matched(i) = true;
    end
end
initiallyMatched = sum(matched);
fprintf(' → Initially matched diet metabolites: %d / %d\n', initiallyMatched, totalDietRxns);

% prepare tissue lists
tissuePrefixes = {'AD_', 'SK_', 'GN_', 'OO_', 'EN_'};
bloodConnectedTissues = {'AD_', 'SK_', 'GN_', 'EN_'};

% handle unmatched diet metabolites
unmatchedIdx = find(~matched);
metsNotInModel = {};
fprintf('\n→ Handling unmatched metabolites (%d total)...\n', numel(unmatchedIdx));
for k = 1:numel(unmatchedIdx)
    i = unmatchedIdx(k);
    dietMet = dietMets{i};
    bloodMet = dietMet;
    rxnDietName = ['Diet_EX_' dietMet];
    % rxnExitName = ['Exit_EX_' dietMet];

    % print progress for each unmatched metabolite
    fprintf('   [%d/%d] Processing %s ... ', k, numel(unmatchedIdx), dietMet);

    % if blood metabolite already in model
    if ismember(bloodMet, model.mets)
        fprintf('found in blood.\n');
        if ~ismember(rxnDietName, model.rxns)
            model = addExchangeRxn(model, bloodMet, -(fluxValues(i)/(24*10000)), 0);
            model.rxns{end} = rxnDietName;
            % model.ub(end) = 0;
        else
            idxRxn = find(strcmp(model.rxns, rxnDietName),1);
            model.lb(idxRxn) = -(fluxValues(i)/(24*10000));
            % model.ub(idxRxn) = 0;
        end
        % if ~ismember(rxnExitName, model.rxns)
        %     model = addReaction(model, rxnExitName, 'metaboliteList', {bloodMet}, 'stoichCoeffList', -1, 'lowerBound', 0, 'upperBound', 1000);
        % end
        continue;
    end

    foundInTissue = false;
    foundInOocyte = false;
    modelMetsClean = regexprep(model.mets, '\[[a-zA-Z]\]', '');
    matchedTissues = {};
    matchedTissuePrefixes = {};

    % ===== 1. FIND MATCHES IN ALL TISSUES =====
    for t = 1:numel(tissuePrefixes)
        tissuePrefix = tissuePrefixes{t};
        dietMetBase = regexprep(dietMet, '\[[a-zA-Z]\]', '');

        % [e] version
        tissueMetE = [tissuePrefix dietMetBase '[e]'];
        if ismember(tissueMetE, model.mets)
            tissueMet = tissueMetE;
            if strcmp(tissuePrefix, 'OO_')
                foundInOocyte = true;
            else
                foundInTissue = true;
            end
            matchedTissues{end+1} = tissueMet;
            matchedTissuePrefixes{end+1} = tissuePrefix;
            continue;
        end

        % non-[e] version
        tissueMetBaseClean = [tissuePrefix dietMetBase];
        nonE_idx = find(startsWith(modelMetsClean, tissueMetBaseClean));
        nonE_idx = nonE_idx(~endsWith(model.mets(nonE_idx), '[e]'));
        if ~isempty(nonE_idx)
            tissueMet = model.mets{nonE_idx(1)};
            if strcmp(tissuePrefix, 'OO_')
                foundInOocyte = true;
            else
                foundInTissue = true;
            end
            matchedTissues{end+1} = tissueMet;
            matchedTissuePrefixes{end+1} = tissuePrefix;
        end
    end

    % ===== 2. ADD REACTIONS (ONLY ONCE) =====
    if foundInTissue || foundInOocyte
        if ~ismember(bloodMet, model.mets)
            model = addMetabolite(model, bloodMet);
        end

        addedTransport = 0;
        addedExchange = 0;

        for m = 1:numel(matchedTissues)
            tissueMet = matchedTissues{m};
            tissuePrefix = matchedTissuePrefixes{m};
            dietMetBase = regexprep(dietMet, '\[[a-zA-Z]\]', '');

            if ismember(tissuePrefix, bloodConnectedTissues)
                transportRxn = ['Tr' dietMetBase '_' tissuePrefix(1:end-1) '2Bl'];
                if ~ismember(transportRxn, model.rxns)
                    model = addReaction(model, transportRxn, ...
                        'metaboliteList', {bloodMet, tissueMet}, ...
                        'stoichCoeffList', [-1, 1], ...
                        'lowerBound', -1000000, 'upperBound', 1000000);
                    addedTransport = addedTransport + 1;
                end

            elseif strcmp(tissuePrefix, 'OO_')
                gnMetE = ['GN_' dietMetBase '[e]'];
                gnMetC = ['GN_' dietMetBase '[c]'];

                if ismember(gnMetE, model.mets)
                    gnMet = gnMetE;
                elseif ismember(gnMetC, model.mets)
                    gnMet = gnMetC;
                else
                    fprintf('   GN_ metabolite for %s not found.\n', dietMetBase);
                    continue;
                end

                trOO2GN = ['Tr' dietMetBase '_OO2GN'];
                if ~ismember(trOO2GN, model.rxns)
                    model = addReaction(model, trOO2GN, ...
                        'metaboliteList', {gnMet, tissueMet}, ...
                        'stoichCoeffList', [-1, 1], ...
                        'lowerBound', -1000000, 'upperBound', 1000000);
                    addedTransport = addedTransport + 1;
                end

                trGN2Bl = ['Tr' dietMetBase '_GN2Bl'];
                if ~ismember(trGN2Bl, model.rxns)
                    model = addReaction(model, trGN2Bl, ...
                        'metaboliteList', {bloodMet, gnMet}, ...
                        'stoichCoeffList', [-1, 1], ...
                        'lowerBound', -1000000, 'upperBound', 1000000);
                    addedTransport = addedTransport + 1;
                end
            end
        end

        % Exchange reactions
        if ~ismember(rxnDietName, model.rxns)
            model = addExchangeRxn(model, bloodMet, -(fluxValues(i)/(24*10000)), 0);
            model.rxns{end} = rxnDietName;
            % model.ub(end) = 0;
            addedExchange = addedExchange + 1;
        else
            idxRxn = find(strcmp(model.rxns, rxnDietName), 1);
            model.lb(idxRxn) = -(fluxValues(i)/(24*10000));
            % model.ub(idxRxn) = 0;
        end

        % if ~ismember(rxnExitName, model.rxns)
        %     model = addReaction(model, rxnExitName, ...
        %         'metaboliteList', {bloodMet}, ...
        %         'stoichCoeffList', -1, ...
        %         'lowerBound', 0, 'upperBound', 1000);
        %     addedExchange = addedExchange + 1;
        % end

        fprintf('   Added %d transport and %d exchange reactions for %s.\n', ...
            addedTransport, addedExchange, dietMet);
    else
        fprintf('not found in any tissue or blood.\n');
        metsNotInModel{end+1} = dietMet;
    end
end


% unmatched count
unmatchedCount = numel(metsNotInModel);

fprintf('\n========== MATCHING SUMMARY ==========\n');
fprintf(' → Initially matched: %d\n', initiallyMatched);
fprintf(' → Matched after tissue/blood checks: %d\n', totalDietRxns - unmatchedCount);
fprintf(' → Unmatched metabolites remaining: %d\n', unmatchedCount);
if unmatchedCount > 0
    fprintf(' → List of unmatched metabolites:\n');
    disp(metsNotInModel');
end
fprintf('=====================================\n\n');
end
