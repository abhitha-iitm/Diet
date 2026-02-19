function [newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges] = nutritionAlgorithm_new(model,rois,options)
% This algorithm identifies the minimal changes to a diet necessary to get
% a desired change in one or more reactions of interest (rois). If a
% metabolite is entered instead of a reaction, the algorithm will optimize
% the diet with a sink or demand reaction for the corresponding metabolite
% of interest. For a walkthrough of the algorithm, see NutritionAlgorithmWalkthrough.mlx
% To cite the algorithm, please cite Weston and Thiele, 2022 and the COBRA
% Toolbox as specified on opencobra.github.io/cobratoolbox/stable/cite.html
% USAGE:
%
%    [newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges,detailedAnalysis] = nutritionAlgorithm(model,rois,roisMinMax,options)
%
% INPUTS:
%    model:          COBRA model structure with minimal fields:
%                      * .S
%                      * .c
%                      * .ub
%                      * .lb
%                      * .mets
%                      * .rxns
%   rois:          cell array of all reactions of interest
%   roisMinMax:    cell array of 'min'/'max' entries for rois
%
% OPTIONAL INPUTS:
%   options:  Structure containing the optional specifications:
%
%       * .display: display results "off" or "on"?
%
%       * .roiWeights:   a vector of weights for each reaction of interest
%       default is equal to 1
%
%       * .targetedDietRxns: A nx2 cell array that specifies any dietary
%       items to target and the corresponding weight for adding the item.
%
%       * .foodRemovalWeighting: Determines the relationship of food
%       removal weight with the weights specified by targetedDietRxns.
%       The following are valid inputs for foodRemovalWeighting.
%            - 'ones'    -> all dietary reactions from targetedDietRxns have a removal weight of one
%            - 'inverse' -> all dietary reactions from targetedDietRxns have a removal weight of one
%            - 'ditto'   -> weights are equal to that of targetedDietRxns
%            added weights
%            - nx2 cell array -> a customized cell array that functions like
%            targedDietRxns but instead allows costomized weights for
%            removing a food item rather than adding it.
%            - {} empty cell array -> (default) An empty array
%            assumes all dietary reactions are available from removal with
%            a weight of one
%
%       * .slnType: Specify if solution should be 'Detailed' or 'Quick'.
%                   Default setting is 'Detailed'
%
%       * .roiBound: 'Unbounded' or 'Bounded'. Default is 'Bounded'.
%
%       * .foodAddedLimit: Specify a limit to the points produced by
%                     adding food to the diet
%
%       * .foodRemovedLimit: Specify a limit to the points produced by
%                     removing food from the diet
%
%       * .OFS: the Objective Flux Scalar initiates a limiting threshold
%       for the solution's objective function performance. A OFS of 1 means
%       that the nutrition algorithm solution will produce a result that is
%       atleast equal to, or greater than the maximum flux of the objective
%       reaction on the original diet
%
%
% OUTPUT:
%    newDietModel:   An copy of the input model with updated diet
%                    reaction bounds to reflect recomended dietary changes
%
%   pointsModel:     The resulting model that is used to identify
%                    recomended dietary changes. It includes points
%                    reactions and food added/removed reactions.
%
%   roiFlux:         Returns the flux values for each roi in the points sln
%
%   pointsModelSln:  Returns the entire points solution to pointsModel
%
%   menuChanges:     Summarizes the recommended dietary changes
%
%   detailedAnalysis: Provides solutions for each simulation conducted in
%                     the detailed analysis
%
% .. Authors: - Bronson R. Weston   2022

% Determine if any rois are metabolites
obj=model.rxns(find(model.c==1));
if isfield(model, 'osenseStr')
    objMinMax=model.osenseStr;
    if strcmp(objMinMax,'min')
        error('osensStr is set to "min". A minimized objective function is not currently supported by the nutrition algorithm. It is recomended that you flip the reaction such that products<-reactants and then set osensStr to maximize')
    end
else
    error('Model fieldname "osenseStr" required to specify if the objective function is to be minimized or maximized.')
end
if isa(rois,'char')
    rois={rois};
end


%If any reactions are targeting the exit reactions accumulation of a specific element
%then create relative changes to the model
for i=1:length(rois)
    if contains(rois{i},'Exit(')

    end
end

% Added targetFlux optional input
targetFlux = [];
if exist('options','var') && isfield(options,'targetFlux')
    targetFlux = options.targetFlux;
    if length(targetFlux) ~= length(rois)
        error('Length of targetFlux must match the number of ROIs');
    end
end

objIndex = find(model.c ~= 0);
if isempty(objIndex)
    error('No objective function set in model.');
end
obj = model.rxns{objIndex(1)};

if isfield(model, 'osenseStr')
    objMinMax = model.osenseStr;
    if strcmp(objMinMax,'min')
        error('osensStr is set to "min". Change to maximize.');
    end
else
    error('Model fieldname "osenseStr" required.');
end


%initialize optional variables
roiWeights=ones(1,length(rois));
targetedDietRxns={};
slnType='Detailed';
roiBound='Bounded';
foodAddedLimit=1000;
foodRemovedLimit=1000;
foodRemovalWeighting={};
display='on';
OFS=1;
if exist('options','var') && ~isempty(options)
    fn = fieldnames(options);
    for k=1:numel(fn)
        if strcmp(fn{k},'roiWeights')
            roiWeights=options.roiWeights;
            if length(roiWeights)~=length(rois)
                error('The length of the roiWeights vector must be equivalent to the length of rois.')
            end
            if any(roiWeights<0)
                error('All roiWeight elements must be greater than zero.')
            end
        elseif strcmp(fn{k},'OFS')
            if OFS>1 || OFS<0
                OFS=1;
                warning('Invalid OFS specified. OFS set to 1.')
            else
                OFS=options.OFS;
            end
        elseif strcmp(fn{k},'foodRemovalWeighting')
            foodRemovalWeighting=options.foodRemovalWeighting;

            if ischar(foodRemovalWeighting) && ~strcmp(foodRemovalWeighting,'ones') && ~strcmp(foodRemovalWeighting,'inverse') ...
                    && ~strcmp(foodRemovalWeighting,'ditto') %&& ~strcmp(foodRemovalWeighting,'cost')
                error('Invalid foodRemovalWeighting input. Must be "ones", "inverse", "ditto" or a nx2 cell array specifiying specific deitary reactions (first column) and the associated weight (second column)')
            end
        elseif strcmp(fn{k},'roiBound')
            roiBound=options.roiBound;
            if ~strcmp(roiBound,'Unbounded') && ~strcmp(roiBound,'Bounded')
                error('Invalid roiBound input. Must be "Unbounded" or "Bounded"')
            end
        elseif strcmp(fn{k},'display')
            display=options.display;
            if ~strcmp(display,'on') && ~strcmp(display,'off')
                error('Invalid display input. Must be "on" or "off"')
            end
        elseif strcmp(fn{k},'foodAddedLimit')
            foodAddedLimit=options.foodAddedLimit;
        elseif strcmp(fn{k},'foodRemovedLimit')
            foodRemovedLimit=options.foodRemovedLimit;
            if ~isnumeric(foodRemovedLimit) || foodRemovedLimit<0
                error('Invalid foodRemovedLimit input')
            end
        elseif strcmp(fn{k},'slnType')
            slnType=options.slnType;
            if strcmp(slnType,'Quick')
                detailedAnalysis=[];
            end
            if ~strcmp(slnType,'Detailed') && ~strcmp(slnType,'Quick')
                error('Invalid slnType input. Must be "Detailed" or "Quick"')
            end
        elseif strcmp(fn{k},'targetedDietRxns')
            targetedDietRxns=options.targetedDietRxns;
        elseif strcmp(fn{k},'targetFlux')
            targetFlux = options.targetFlux;
            if length(targetFlux) ~= length(rois)
                error('Length of targetFlux must match the number of ROIs');
            end
        elseif strcmp(fn{k},'pcosFlux')
            pcosFlux = options.pcosFlux;
            if length(pcosFlux) ~= length(rois)
                error('Length of pcosFlux must match the number of ROIs');
            end
        else
            error(['Invalid "options" field entered: ', fn{k}])
        end
    end
end


if strcmp(display,'on')
    disp('_____________________________________________________')
end

newDietModel=model; %Copy original instance of model for new diet pointsModel
pointsModel=model; %Copy original instance of model for points pointsModel


%Calculate newDietModel objective function and stores the indices of roi reactions
objIndex=find(model.c==1);

roiIndexO=zeros(1,length(rois));
for i=1:length(roiIndexO)
    roiInd = find(strcmp(newDietModel.rxns,rois{i}));
    if isempty(roiInd)
        error(['The following roi is not a valid rxn in the model: ', rois{i}])
    end
    roiIndexO(i)=find(strcmp(newDietModel.rxns,rois{i}));
end


pointsModel=addMetabolite(pointsModel, 'unitOfFoodAdded[dP]');
pointsModel=addMetabolite(pointsModel, 'unitOfFoodRemoved[dP]');
pointsModel=addMetabolite(pointsModel, 'unitOfFoodChange[dP]');
pointsModel=addMetabolite(pointsModel, 'roiPoint[roiP]');
pointsModel=addMetabolite(pointsModel, 'point[P]');

nPrint = min(10, length(rois));
for i = 1:nPrint

    roiName = char(rois{i});
    tgt     = targetFlux(i);
    pcos    = options.pcosFlux(i);
    delta   = pcos - tgt;

    fprintf('\n===== DEBUG ROI %d =====\n', i);
    fprintf('ROI Name  : %s\n', roiName);
    fprintf('Target    : %.4f\n', tgt);
    fprintf('PCOS      : %.4f\n', pcos);
    fprintf('Delta     : %.4f\n', delta);

end

if ~isempty(targetFlux)

    newRois    = {};
    newWeights = [];

    for i = 1:length(rois)

        roiName = char(rois{i});
        tgt     = targetFlux(i);
        pcos    = options.pcosFlux(i);
        delta   = pcos - tgt;

        % POINTS MODEL DECOMPOSITION

        rxnIdxP = find(strcmp(pointsModel.rxns, roiName), 1);
        if isempty(rxnIdxP)
            error(['ROI reaction not found: ', roiName]);
        end

        stoichP   = full(pointsModel.S(:, rxnIdxP));
        metsMaskP = (stoichP ~= 0);
        metsListP = pointsModel.mets(metsMaskP)';
        coeffP    = stoichP(metsMaskP);

        % Remove original reaction
        evalc('[pointsModel,~,~] = removeRxns(pointsModel, {roiName});');

        baseP = [roiName '_base'];
        dposP = [roiName '_dev_pos'];
        dnegP = [roiName '_dev_neg'];

        SnewP = zeros(numel(metsListP),3);
        SnewP(:,1) = coeffP;
        SnewP(:,2) = coeffP;
        SnewP(:,3) = -coeffP;

        % Deviation bounds
        if delta > 0
            devPosLB = 0;
            devPosUB = delta;
            devNegLB = 0;
            devNegUB = 0;
        elseif delta < 0
            devPosLB = 0;
            devPosUB = 0;
            devNegLB = 0;
            devNegUB = abs(delta);
        else
            devPosLB = 0;
            devPosUB = 0;
            devNegLB = 0;
            devNegUB = 0;
        end

        pointsModel = addMultipleReactions(pointsModel, ...
            {baseP, dposP, dnegP}, ...
            metsListP, ...
            SnewP, ...
            'lb', [tgt, devPosLB, devNegLB], ...
            'ub', [tgt, devPosUB, devNegUB]);


        % base + dev_pos - dev_neg = PCOS flux

        linkMet = [roiName '_pcos_balance'];

        pointsModel = addMetabolite(pointsModel, linkMet);

        linkIdx = strcmp(pointsModel.mets, linkMet);

        % base produces link metabolite
        pointsModel.S(linkIdx, strcmp(pointsModel.rxns, baseP)) = 1;

        % dev_pos produces link metabolite
        pointsModel.S(linkIdx, strcmp(pointsModel.rxns, dposP)) = 1;

        % dev_neg consumes link metabolite
        pointsModel.S(linkIdx, strcmp(pointsModel.rxns, dnegP)) = -1;

        % fixed PCOS demand
        pointsModel = addReaction(pointsModel, ...
            [roiName '_PCOS_total'], ...
            {linkMet}, ...
            -1, ...
            false, ...
            pcos, ...
            pcos);

        % Store dev reactions for ROI penalty
        newRois{end+1} = dposP;
        newRois{end+1} = dnegP;

        if exist('roiWeights','var') && length(roiWeights) >= i
            newWeights = [newWeights, roiWeights(i), roiWeights(i)];
        else
            newWeights = [newWeights, 1, 1];
        end

    end

    rois       = newRois;
    roiWeights = newWeights;

    % ROI PENALTY SYSTEM

    roiPointIdx = find(strcmp(pointsModel.mets,'roiPoint[roiP]'));

    for i = 1:length(rois)
        devIdx = find(strcmp(pointsModel.rxns, rois{i}));
        if ~isempty(devIdx)
            pointsModel.S(roiPointIdx, devIdx) = roiWeights(i);
        end
    end

    % Pool ROI points into global point[P]
    pointsModel = addReaction(pointsModel,...
        'Point_EX_roiPoints[roiP]_[P]',...
        {'roiPoint[roiP]','point[P]'},...
        [-1 1],...
        true,...
        -1000000,...
        1000000);

end


%add all diet exchange reactions to targetedDietRxns
if isempty(targetedDietRxns)
    %Add all Diet_EX reactions, and set weight to 1
    dietRxns=find(contains(pointsModel.rxns,'Diet_EX'));
    targetedDietRxns=[pointsModel.rxns(dietRxns),num2cell(ones(length(dietRxns),1))];
elseif any(strcmp(targetedDietRxns(:,1),'All')) && length(targetedDietRxns(:,1))>1
    dietRxns=find(contains(pointsModel.rxns,'Diet_EX'));
    targetedFoodItemsTemp=[pointsModel.rxns(dietRxns), ...
        num2cell(cell2mat(targetedDietRxns(strcmp(targetedDietRxns(:,1),'All'),2))*ones(length(dietRxns),1))];
    [~,ai,bi]=intersect(targetedFoodItemsTemp(:,1),targetedDietRxns(:,1));
    targetedFoodItemsTemp(ai,2)=targetedDietRxns(bi,2);
    targetedDietRxns=targetedFoodItemsTemp;
elseif any(strcmp(targetedDietRxns(:,1),'All'))
    dietRxns=find(contains(pointsModel.rxns,'Diet_EX'));
    targetedDietRxns=[pointsModel.rxns(dietRxns),num2cell(targetedDietRxns{1,2}*ones(length(dietRxns),1))];
end

%Set up foodRemovalWeighting
if isempty(foodRemovalWeighting)
    %Add all Diet_EX reactions, and set weight to 1
    dietRxns=find(contains(pointsModel.rxns,'Diet_EX'));
    foodRemovalWeighting=[pointsModel.rxns(dietRxns),num2cell(ones(length(dietRxns),1))];
elseif ischar(foodRemovalWeighting)
    foodRemoveTemp=targetedDietRxns;
    switch foodRemovalWeighting
        case 'ones'
            foodRemoveTemp(:,2)=num2cell(ones(length(foodRemoveTemp(:,2)),1));
        case 'ditto'
        case 'inverse'
            foodRemoveTemp(:,2)=num2cell(cell2mat(targetedDietRxns(:,2)).^-1);
            %         case 'cost'
            %             foodRemoveTemp(:,2)=num2cell(cell2mat(targetedDietRxns(:,2))*-1);
    end
    foodRemovalWeighting=foodRemoveTemp; clear foodRemoveTemp;
elseif any(strcmp(foodRemovalWeighting(:,1),'All')) && length(foodRemovalWeighting(:,1))>1
    dietRxns=find(contains(pointsModel.rxns,'Diet_EX'));
    targetedFoodItemsTemp=[pointsModel.rxns(dietRxns), ...
        num2cell(cell2mat(foodRemovalWeighting(strcmp(foodRemovalWeighting(:,1),'All'),2))*ones(length(dietRxns),1))];
    [~,ai,bi]=intersect(targetedFoodItemsTemp(:,1),foodRemovalWeighting(:,1));
    targetedFoodItemsTemp(ai,2)=foodRemovalWeighting(bi,2);
    foodRemovalWeighting=targetedFoodItemsTemp;
elseif any(strcmp(foodRemovalWeighting(:,1),'All'))
    dietRxns=find(contains(pointsModel.rxns,'Diet_EX'));
    foodRemovalWeighting=[pointsModel.rxns(dietRxns),num2cell(foodRemovalWeighting{1,2}*ones(length(dietRxns),1))];
end


%Add "Food_Added" reactions to points model
Mets= pointsModel.mets;
[~,ai,bi]=intersect(pointsModel.rxns,targetedDietRxns(:,1));
if isempty(bi)
    error('targetedDietRxns does not include any valid reactions in the model')
end
foodRxns=pointsModel.rxns(ai);
foodRxns=regexprep(foodRxns,'Diet_EX_','Food_Added_EX_');
sMatrix=pointsModel.S(:,ai);
f=find(strcmp(pointsModel.mets,'unitOfFoodAdded[dP]'));
sMatrix(f,:)=-1*cell2mat(targetedDietRxns(bi,2)).';
pointsModel = addMultipleReactions(pointsModel, foodRxns, Mets, sMatrix, 'lb', -1000000*ones(1,length(foodRxns)), 'ub', zeros(1,length(foodRxns)));
pointsModel = addMultipleReactions(pointsModel, {'Point_EX_unitOfFoodRemoved2Change[dp]','Point_EX_unitOfFoodAdded2Change[dp]','Point_EX_unitOfFoodChange[dP]_[P]','Point_EX_Point[P]'}, {'unitOfFoodRemoved[dP]','unitOfFoodAdded[dP]','unitOfFoodChange[dP]','point[P]'}, [-1 0 0 0;0 -1 0 0;1 1 -1 0;0 0 1 -1], 'lb', [-1000000,-1000000, -1000000,-1000000], 'ub', [foodRemovedLimit,foodAddedLimit,1000000,1000000]);


%Add "Food Removed" reactions to points model
[~,ai,bi]=intersect(pointsModel.rxns,foodRemovalWeighting(:,1));
bi=bi(pointsModel.lb(ai)<0); % only includes removal reactions for dietary reactions that have a non-zero influx
ai=ai(pointsModel.lb(ai)<0);
if ~isempty(ai)
    foodRxns=pointsModel.rxns(ai);
    foodRxns=regexprep(foodRxns,'Diet_EX_','Food_Removed_EX_');
    sMatrix=-1*pointsModel.S(:,ai);
    f=find(strcmp(pointsModel.mets,'unitOfFoodRemoved[dP]'));
    sMatrix(f,:)=-1*cell2mat(foodRemovalWeighting(bi,2)).';
    pointsModel = addMultipleReactions(pointsModel, foodRxns, pointsModel.mets, sMatrix, 'lb', pointsModel.lb(ai), 'ub', zeros(1,length(ai)));
end

disp('Checking Targeted Diet Reactions and Food Added/Removed Reactions:');

% Display dietary reactions and weights
for i = 1:size(targetedDietRxns,1)
    disp([targetedDietRxns{i,1}, ' weight: ', num2str(targetedDietRxns{i,2})]);
end

% Check food added reactions in pointsModel
foodAddedRxns = regexprep(targetedDietRxns(:,1), 'Diet_EX_', 'Food_Added_EX_');
foodAddedPresent = ismember(foodAddedRxns, pointsModel.rxns);
fprintf('Food_Added reactions present: %d/%d\n', sum(foodAddedPresent), length(foodAddedRxns));

% Check food removed reactions in pointsModel
foodRemovedRxns = regexprep(targetedDietRxns(:,1), 'Diet_EX_', 'Food_Removed_EX_');
foodRemovedPresent = ismember(foodRemovedRxns, pointsModel.rxns);
fprintf('Food_Removed reactions present: %d/%d\n', sum(foodRemovedPresent), length(foodRemovedRxns));


%Find solution
devIdx  = find(strcmp(pointsModel.rxns,'Point_EX_roiPoints[roiP]_[P]'));
foodIdx = find(strcmp(pointsModel.rxns,'Point_EX_unitOfFoodChange[dP]_[P]'));

Wdev  = 1;
Wfood = 0;

pointsModel.c(:) = 0;
pointsModel.c(devIdx)  = Wdev;
pointsModel.c(foodIdx) = Wfood;
pointsModel.osenseStr = 'min';

pointsModelSln = optimizeCbModel(pointsModel);

fprintf('\n================ DEBUG START ================\n');

fprintf('\n--- REACTION EXISTENCE CHECK ---\n');

baseIdx = find(strcmp(pointsModel.rxns, baseP));
dposIdx = find(strcmp(pointsModel.rxns, dposP));
dnegIdx = find(strcmp(pointsModel.rxns, dnegP));

fprintf('Base reaction index : %s\n', mat2str(baseIdx));
fprintf('Dev+ reaction index : %s\n', mat2str(dposIdx));
fprintf('Dev- reaction index : %s\n', mat2str(dnegIdx));

fprintf('\n--- DUPLICATE CHECK ---\n');
fprintf('Number of base reactions : %d\n', length(baseIdx));
fprintf('Number of dev+ reactions : %d\n', length(dposIdx));
fprintf('Number of dev- reactions : %d\n', length(dnegIdx));

fprintf('\n--- BOUNDS CHECK ---\n');
if ~isempty(baseIdx)
    fprintf('Base LB : %s\n', mat2str(pointsModel.lb(baseIdx)));
    fprintf('Base UB : %s\n', mat2str(pointsModel.ub(baseIdx)));
end

if ~isempty(dposIdx)
    fprintf('Dev+ LB : %s\n', mat2str(pointsModel.lb(dposIdx)));
    fprintf('Dev+ UB : %s\n', mat2str(pointsModel.ub(dposIdx)));
end

if ~isempty(dnegIdx)
    fprintf('Dev- LB : %s\n', mat2str(pointsModel.lb(dnegIdx)));
    fprintf('Dev- UB : %s\n', mat2str(pointsModel.ub(dnegIdx)));
end

fprintf('\n--- STOICHIOMETRY CHECK ---\n');
if ~isempty(baseIdx)
    sto = full(pointsModel.S(:, baseIdx(1)));
    fprintf('Base stoich nonzero count : %d\n', sum(sto~=0));
end

if ~isempty(dposIdx)
    sto = full(pointsModel.S(:, dposIdx(1)));
    fprintf('Dev+ stoich nonzero count : %d\n', sum(sto~=0));
end

if ~isempty(dnegIdx)
    sto = full(pointsModel.S(:, dnegIdx(1)));
    fprintf('Dev- stoich nonzero count : %d\n', sum(sto~=0));
end

fprintf('\n--- QUICK SOLVER TEST ---\n');
testSol = optimizeCbModel(pointsModel);
disp('--- CHECK DIET BOUNDS ---')
rxn = 'Diet_EX_arach[e]';

lb = pointsModel.lb(strcmp(pointsModel.rxns, rxn));
ub = pointsModel.ub(strcmp(pointsModel.rxns, rxn));

fprintf('%s LB = %g\n', rxn, lb);
fprintf('%s UB = %g\n', rxn, ub);


if isempty(testSol.f) || isnan(testSol.f)
    fprintf('Solver returned infeasible solution.\n');
else
    fprintf('Objective value : %.6f\n', testSol.f);

    if ~isempty(baseIdx)
        fprintf('Base flux : %s\n', mat2str(testSol.v(baseIdx)));
    end
    if ~isempty(dposIdx)
        fprintf('Dev+ flux : %s\n', mat2str(testSol.v(dposIdx)));
    end
    if ~isempty(dnegIdx)
        fprintf('Dev- flux : %s\n', mat2str(testSol.v(dnegIdx)));
    end
end

fprintf('\n--- ROI POINTS CHECK ---\n');
roiPointIdx = find(strcmp(pointsModel.rxns,'Point_EX_roiPoints[roiP]_[P]'));
if ~isempty(roiPointIdx)
    fprintf('ROI Points flux : %.6f\n', testSol.v(roiPointIdx));
else
    fprintf('ROI Points reaction not found.\n');
end

fprintf('\n================ DEBUG END ==================\n\n');


dietFlux = pointsModelSln.v(foodIdx);
devFlux  = pointsModelSln.v(devIdx);

% Compute weighted score ONLY for reporting
weightedTotalPoints = (Wdev * devFlux) + (Wfood * dietFlux);

if strcmp(display,'on')
    disp(['Weighted Total Points = ', num2str(weightedTotalPoints)]);
    disp(['Solution points (LP objective) = ', num2str(pointsModelSln.f)]);
    disp([num2str(dietFlux), ' come from diet']);
    disp([num2str(devFlux),  ' come from rois']);
end

foodAddedIndexes=find(contains(pointsModel.rxns,'Food_Added_EX_'));
foodRemovedIndexes=find(contains(pointsModel.rxns,'Food_Removed_EX_'));
slnIndexes1=foodAddedIndexes(pointsModelSln.v(foodAddedIndexes)<0);
slnIndexes2=foodRemovedIndexes(pointsModelSln.v(foodRemovedIndexes)<0);

menuChanges=table([pointsModel.rxns(slnIndexes1);pointsModel.rxns(slnIndexes2)],pointsModelSln.v([slnIndexes1;slnIndexes2]),'VariableNames',{'Food Rxn', 'Flux'});
if strcmp(display,'on')
    menuChanges
end

%Add and remove relevant food items from diet in newDietModel
foodItemsAdd= regexprep(pointsModel.rxns(slnIndexes1),'Food_Added_EX_','Diet_EX_');
foodItemsRemove= regexprep(pointsModel.rxns(slnIndexes2),'Food_Removed_EX_','Diet_EX_');
modelOindexAdd=zeros(1,length(foodItemsAdd));
sl2IndexAdd=zeros(1,length(foodItemsAdd));
modelOindexRemove=zeros(1,length(foodItemsRemove));
sl2IndexRemove=zeros(1,length(foodItemsRemove));
for i=1:length(foodItemsAdd)
    modelOindexAdd(i)=find(strcmp(newDietModel.rxns,foodItemsAdd(i)));
    sl2IndexAdd(i)=find(strcmp(pointsModel.rxns,foodItemsAdd(i)));
end
for i=1:length(foodItemsRemove)
    modelOindexRemove(i)=find(strcmp(newDietModel.rxns,foodItemsRemove(i)));
    sl2IndexRemove(i)=find(strcmp(pointsModel.rxns,foodItemsRemove(i)));
end
% newDietModel.lb(modelOindexAdd)=(pointsModelSln.v(sl2IndexAdd)+pointsModelSln.v(slnIndexes1))*1.01;
newDietModel.lb(modelOindexAdd)=(pointsModelSln.v(sl2IndexAdd)+pointsModelSln.v(slnIndexes1));
newDietModel.ub(modelOindexAdd)=(pointsModelSln.v(sl2IndexAdd)+pointsModelSln.v(slnIndexes1));
% newDietModel.lb(modelOindexRemove)=(pointsModelSln.v(sl2IndexRemove)-pointsModelSln.v(slnIndexes2))*1.01;
newDietModel.lb(modelOindexRemove)=(pointsModelSln.v(sl2IndexRemove)-pointsModelSln.v(slnIndexes2));
newDietModel.ub(modelOindexRemove)=(pointsModelSln.v(sl2IndexRemove)-pointsModelSln.v(slnIndexes2));

if strcmp(display,'on')
    disp('Points Simulation Solution:')
end
for i=1:length(rois)
    ind=find(strcmp(pointsModel.rxns,rois{i}));
    % if strcmp(display,'on')
    %     disp(['   ',rois{i},' flux = ', num2str(pointsModelSln.v(ind(1)))])
    % end
    roiFlux(i)=pointsModelSln.v(ind(1));
end
if strcmp(slnType,'Quick')
    detailedAnalysis=[];
    return
end


%Find new obj flux with new diet
if strcmp(display,'on')
    disp('Detailed Analysis:')
end

% Compute original objective value
sol0 = optimizeCbModel(model);
if isempty(sol0.f) || isnan(sol0.f)
    error('Original model optimization failed');
end
f1 = sol0.f;

% Compute objective with new diet
if model.ub(objIndex) ~= model.lb(objIndex)

    model_Obj = optimizeCbModel(newDietModel);

    if isempty(model_Obj.f) || isnan(model_Obj.f)
        error('New diet model optimization failed');
    end

    f2 = model_Obj.f;

    % Constrain objective to maintain performance
    newDietModel = changeRxnBounds(newDietModel, obj, f2, 'l');  
    % use 'b' instead of 'l' if you want to FIX exactly

else
    f2 = model.ub(objIndex);
end

% Print comparison (only if objective not ROI)
if ~any(strcmp(obj,rois)==1) && strcmp(display,'on')
    disp(['Original objective max flux = ', num2str(f1), ...
          '  |  New objective max flux = ', num2str(f2)])
end

% Restore original objective bounds & sense
newDietModel.ub(objIndex) = model.ub(objIndex);
newDietModel.lb(objIndex) = model.lb(objIndex);
newDietModel.osenseStr = objMinMax;

if strcmp(display,'on')
    disp('_____________________________________________________')
end
detailedAnalysis = [];