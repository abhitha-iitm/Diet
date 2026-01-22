function [model] = convert_EX_to_diet(model)
%
% convert_EX_to_diet takes a model with typical exchange reactions 
% (i.e., 'EX_') and translates them to uptake (diet) and secretion (exit) 
% reactions.
%
% USAGE:
%
% [model] = convert_EX_to_diet(model)
%
% INPUTS:
%    model:          COBRA model structure with minimal fields:
%                      * .S
%                      * .c
%                      * .ub
%                      * .lb
%                      * .mets  
%                      * .rxns   
%
%Outputs
%   model: Returns the input model with new uptake ('Diet_') and secretion 
%          ('Exit_') reactions
%
%Authors: Bronson R. Weston 2022


exRxns=[];
for i=1:length(model.rxns)
    tmp=strsplit(model.rxns{i},'EX_');
    if isempty(tmp{1})
        exRxns=[exRxns,i];
    end 
end

% Use S matrix to determine which metabolites are in exchange reactions
S=model.S;
exMets=[];
for i=1:length(exRxns)
    x=find(S(:,exRxns(i))~=0).';
    exMets=[exMets,x];

end
exMets=unique(exMets);

% Create environmental metabolites for each exchange metabolite
if any(contains(model.rxns,'Diet_')) 
    error('Model already contains dietary reactions')
end
envMets=model.mets(exMets);

% Create new exchange reactions importing metabolites in dietary reactions.
% Set ub and lb according to lb of original exchange reaction
%ubs=model.ub(exRxns);
%lbs=model.lb(exRxns);
model.rxns(exRxns)=strcat('Diet_',model.rxns(exRxns));
%model.lb(exRxns)=model.lb(exRxns);
%model.ub(exRxns)=model.lb(exRxns);

% % Add an exit reaction for each original EX_ reaction, matching Diet_EX_ naming style
% for i = 1:length(exRxns)
%     % Get the original EX_ reaction name (without Diet_ prefix)
%     origRxnName = model.rxns{exRxns(i)};
%     if startsWith(origRxnName, 'Diet_')
%         origRxnName = origRxnName(6:end); % remove 'Diet_' prefix
%     end
%     % exitRxnName = ['Exit_' origRxnName];
%     % metIdx = find(S(:,exRxns(i)) ~= 0);
%     % metName = model.mets{metIdx};
%     % model = addMultipleReactions(model, {exitRxnName}, {metName}, -1, 'lb', 0, 'ub', 1000);
% end


% === Count and display Exchange, Diet, and Exit reactions ===

% Find all types of reactions by name
exRxns = find(contains(model.rxns, 'EX_') & ~contains(model.rxns, 'Diet_'));
dietRxns = find(contains(model.rxns, 'Diet_'));
% exitRxns = find(contains(model.rxns, 'Exit_'));

% Display summary counts
fprintf('\n=== Reaction Summary ===\n');
fprintf('Exchange reactions : %d\n', numel(exRxns));
fprintf('Diet reactions     : %d\n', numel(dietRxns));
% fprintf('Exit reactions     : %d\n\n', numel(exitRxns));


