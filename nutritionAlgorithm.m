function [newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges,detailedAnalysis] = nutritionAlgorithm(model,rois,options)


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
if isa(roisMinMax,'char')
    roisMinMax={roisMinMax};
end


%Create a demand or sink reaction, as appropriate, for any rois that are
%metabolites
metRois=[];
for i=1:length(rois)
    if any(strcmp(model.mets,rois{i}))
        metRois=[metRois,i];
        if strcmp(roisMinMax{i},'max')
            model=addDemandReaction(model,rois{i}); %adds demand reaction as 'DM_metabolite'
            rois{i}=['DM_',rois{i}];
            model=changeRxnBounds(model,rois{i},1000000,'u');
        else
            model=addSinkReactions(model,rois(i),-1000000,0);
            rois{i}=['sink_',rois{i}];
        end
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



% Set weighted tissue biomass maintenance objective
tissueBiomassRxns = {'SK_biomass_maintenance','AD_biomass_maintenance','GN_biomass_maintenance','OO_biomass_maintenance','EN_biomass_maintenance'};
weights = [0.2, 0.2, 0.2, 0.2, 0.2];
% tissueBiomassRxns = {'SK_ATPtm','AD_ACCOAC','AD_biomass_maintenance','GN_biomass_maintenance','GN_P450SCC1m','OO_biomass_maintenance','EN_biomass_maintenance','EN_biomass_reaction'};
% weights = [0.2, 0.1,0.1, 0.1, 0.1, 0.2, 0.1, 0.1]; 
model.c(:) = 0;
for i = 1:length(tissueBiomassRxns)
    rxnIdx = find(strcmp(model.rxns, tissueBiomassRxns{i}));
    if ~isempty(rxnIdx)
        model.c(rxnIdx) = weights(i);
    else
        fprintf('Warning: %s not found in model.rxns\n', tissueBiomassRxns{i});
    end
end
objIndex = find(model.c~=0);
if isempty(objIndex) 
    error('No tissue biomass maintenance reactions found in model for objective.');
end
obj = model.rxns{objIndex(1)};
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



% % Ensure model identifiers are cell arrays of char
% if isfield(model,'rxns') && isstring(model.rxns), model.rxns = cellstr(model.rxns); end
% if isfield(model,'mets') && isstring(model.mets), model.mets = cellstr(model.mets); end
% % Determine per-ROI direction without requiring roisMinMax
% % When targetFlux is provided, infer direction from sign of target flux.
% % Otherwise, default to 'min' for all unless options.roiDirection is provided.
% if exist('targetFlux','var') && ~isempty(targetFlux)
%     roiDir = repmat({'min'}, 1, length(rois));
%     for ii=1:length(rois)
%         if targetFlux(ii) > 0
%             roiDir{ii} = 'max';
%         else
%             roiDir{ii} = 'min';
%         end
%     end
% else
%     if exist('options','var') && isfield(options,'roiDirection') && ~isempty(options.roiDirection)
%         roiDir = options.roiDirection;
%         if isstring(roiDir), roiDir = cellstr(roiDir); end
%         if ischar(roiDir), roiDir = repmat({roiDir},1,length(rois)); end
%         if numel(roiDir) ~= length(rois)
%             error('options.roiDirection length must match number of rois');
%         end
%     else
%         roiDir = repmat({'min'}, 1, length(rois));
%     end
% end


% %Create a demand or sink reaction, as appropriate, for any rois that are metabolites
% metRois=[];
% for i=1:length(rois)
%     if any(strcmp(model.mets,rois{i}))
%         metRois=[metRois,i];
%         if strcmp(roiDir{i},'max')
%             model=addDemandReaction(model,rois{i}); %adds demand reaction as 'DM_metabolite'
%             rois{i}=['DM_',rois{i}];
%             model=changeRxnBounds(model,rois{i},1000000,'u');
%         else
%             model=addSinkReactions(model,rois(i),-1000000,0);
%             rois{i}=['sink_',rois{i}];
%         end
%     end
% end


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
        else
            error(['Invalid "options" field entered: ', fn{k}])
        end
    end
end



%adjust ub and lb if roiBound specifies 'Unbound'
% Only adjust bounds when not using targetFlux
if strcmp(roiBound, 'Unbounded') && (isempty(targetFlux))
    for i=1:length(rois)
        f=find(strcmp(model.rxns,rois{i}));
        if strcmp(roiDir{i},'max')
            if model.ub(f)~=0
                model.ub(f)=1000000;
            end
        else
            if model.lb(f)~=0
                model.lb(f)=-1000000;
            end
        end
    end
end

if strcmp(display,'on')
    disp('_____________________________________________________')
end


pointsModel=model; %Copy original instance of model for points pointsModel


tissueBiomassRxns = {'SK_biomass_maintenance','AD_biomass_maintenance','GN_biomass_maintenance','OO_biomass_maintenance','EN_biomass_maintenance'};
weights = [0.2, 0.2, 0.2, 0.2, 0.2]; % adjust as needed
% tissueBiomassRxns = {'SK_ATPtm','AD_ACCOAC','AD_biomass_maintenance','GN_biomass_maintenance','GN_P450SCC1m','OO_biomass_maintenance','EN_biomass_maintenance','EN_biomass_reaction'};
% weights = [0.2, 0.1,0.1, 0.1, 0.1, 0.2, 0.1, 0.1]; 
model.c(:) = 0;
for i = 1:length(tissueBiomassRxns)
    rxnIdx = find(strcmp(model.rxns, tissueBiomassRxns{i}));
    if ~isempty(rxnIdx)  
        model.c(rxnIdx) = weights(i);
    else
        fprintf('Warning: %s not found in model.rxns\n', tissueBiomassRxns{i});
    end
end
sol_CT_weighted = optimizeCbModel(model, 'max');
fprintf('Weighted tissue biomass objective value from  na pc model: %.6f\n', sol_CT_weighted.f);

newDietModel=model; %Copy original instance of model for new diet pointsModel

roiIndexO=zeros(1,length(rois));
for i=1:length(roiIndexO)
    roiInd = find(strcmp(newDietModel.rxns,rois{i}));
    if isempty(roiInd)
        error(['The following roi is not a valid rxn in the model: ', rois{i}])
    end
    roiIndexO(i)=find(strcmp(newDietModel.rxns,rois{i}));
end


% --- Get flux of objective function ---
if ~isempty(objIndex) && any(model.ub(objIndex) ~= model.lb(objIndex))
    % Optimize the new diet model
    model_Obj = optimizeCbModel(newDietModel);
    f1 = model_Obj.f;

    % Print the objective reaction(s) and flux value
    disp(['Initial objective function flux (', strjoin(model.rxns(objIndex), ', '), ') = ', num2str(f1)]);

roiTable = readtable('rois.xlsx');   % <-- your file name
biomassRxns = roiTable{:,1};



% Find their indices in the model
[isPresent, rxnIndices] = ismember(biomassRxns, newDietModel.rxns);


% Get corresponding flux values from the solution
biomassFluxes = model_Obj.x(rxnIndices(isPresent));
disp(biomassFluxes);

% Display results
disp('Flux values for biomass maintenance reactions:');
for i = 1:length(biomassRxns)
    if isPresent(i)
        disp([biomassRxns{i}, ' = ', num2str(model_Obj.x(rxnIndices(i)))]);
    else
        disp([biomassRxns{i}, ' not found in model']);
    end
end


    if ~isnan(f1)
        % --- Step 1: Constrain main objective in pointsModel (optional safeguard)
        pointsModel = changeRxnBounds(pointsModel, model.rxns(objIndex), f1, 'l');

        % --- Step 2: Constrain each tissue biomass reaction to at least control flux ---
        tissueBiomassRxns = {'SK_biomass_maintenance', 'AD_biomass_maintenance', ...
                             'GN_biomass_maintenance', 'OO_biomass_maintenance', ...
                             'EN_biomass_maintenance'};
        % tissueBiomassRxns = {'SK_ATPtm','AD_ACCOAC','AD_biomass_maintenance','GN_biomass_maintenance','GN_P450SCC1m','OO_biomass_maintenance','EN_biomass_maintenance','EN_biomass_reaction'};
        % Use same total flux (not divided)
        target_each = 0.01;

        for i = 1:length(tissueBiomassRxns)
            idx = find(strcmp(pointsModel.rxns, tissueBiomassRxns{i}));
            if ~isempty(idx)
                pointsModel = changeRxnBounds(pointsModel, tissueBiomassRxns{i}, target_each, 'l');
                fprintf('Set lower bound for %s = %.4f\n', tissueBiomassRxns{i}, target_each);
            else
                warning('%s not found in pointsModel.rxns', tissueBiomassRxns{i});
            end
        end
    end

    try
        initRoiFlux = model_Obj.v(roiIndexO);
    catch
        initRoiFlux = NaN(1, length(rois));
    end
else
    f1 = model.ub(objIndex);
    initRoiFlux = NaN(1, length(rois));
end

% --- Display final objective in pointsModel ---
objIdx_points = find(pointsModel.c ~= 0);
if isempty(objIdx_points)
    disp('No objective function is currently set in pointsModel.');
else
    disp('Objective function reaction(s) in pointsModel:');
    disp(pointsModel.rxns(objIdx_points));
    disp('Corresponding coefficients:');
    disp(pointsModel.c(objIdx_points));
end




if isempty(targetFlux)
    % Existing detailed min/max logic remains unchanged
    if strcmp(slnType,'Detailed') && all(~isnan(f1))
        OroiFluxMin=[];
        OroiFluxMax=[];
        for i=1:length(rois)
            if initRoiFlux(i)==newDietModel.lb(roiIndexO(i))
                OroiFluxMin(i)=newDietModel.lb(roiIndexO(i));
                detailedAnalysis.(['Rxn',num2str(i)]).min.OD=model_Obj;
            else
                pointsModel = changeObjective(pointsModel,rois{i});
                pointsModel.osenseStr = 'min';
                if strcmp(rois{i},obj)
                    pointsModel=changeRxnBounds(pointsModel,obj,model.lb(objIndex),'l'); %constrain pointsModel obj flux
                    sln = optimizeCbModel(pointsModel);
                    pointsModel=changeRxnBounds(pointsModel,obj,f1,'l'); %constrain pointsModel obj flux
                else
                    sln = optimizeCbModel(pointsModel);
                end
                detailedAnalysis.(['Rxn',num2str(i)]).min.OD=sln;
                if isnan(sln.f)
                    warning('model input into function does not have a viable initial solution');
                    OroiFluxMin(i)=NaN;
                    OroiFluxMax(i)=NaN;
                    continue
                end
                OroiFluxMin(i)=sln.v(roiIndexO(i));
            end
            if initRoiFlux(i)==newDietModel.ub(roiIndexO(i))
                OroiFluxMax(i)=newDietModel.ub(roiIndexO(i));
                detailedAnalysis.(['Rxn',num2str(i)]).max.OD=model_Obj;
            else
                pointsModel = changeObjective(pointsModel,rois{i});
                pointsModel.osenseStr = 'max';
                if strcmp(rois{i},obj)
                    sln = model_Obj;
                else
                    sln = optimizeCbModel(pointsModel);
                end
                detailedAnalysis.(['Rxn',num2str(i)]).max.OD=sln;
                if isnan(sln.f)
                    warning('Not a viable initial solution, recommend adding nutrients to initial diet');
                    OroiFluxMax(i)=NaN;
                    continue
                end
                OroiFluxMax(i)=sln.v(roiIndexO(i));
            end
        end
    elseif any(isnan(f1))
        for i=1:length(rois)
            OroiFluxMax(i)=NaN;
            OroiFluxMin(i)=NaN;
        end
        warning('Not a viable initial solution, recommend adding nutrients to initial diet');
    end
else
    % When targetFlux is provided, just record current flux
    OroiFluxMin = initRoiFlux;
    OroiFluxMax = initRoiFlux;
    detailedAnalysis = [];
end




% % If targetFlux provided, decompose each ROI into base + deviation reactions
% if ~isempty(targetFlux)
%     newRois = {};
%     newWeights = [];
%     for i = 1:length(rois)
%     roiName = char(rois{i});
%         tgt = targetFlux(i);
%         pcos=biomassFluxes(i);
%         % pointsModel decomposition
%         rxnIdxP = find(strcmp(pointsModel.rxns, roiName), 1);
%         if isempty(rxnIdxP)
%             error(['ROI reaction not found in pointsModel: ', roiName]);
%         end
%         stoichP = full(pointsModel.S(:, rxnIdxP));
%         metsMaskP = (stoichP ~= 0);
%         metsListP = pointsModel.mets(metsMaskP)';
%         coeffP = stoichP(metsMaskP);
%         evalc('[pointsModel,~,~] = removeRxns(pointsModel, {roiName});');
%         baseP = [roiName, '_base'];
%         dposP = [roiName, '_dev_pos'];
%         dnegP = [roiName, '_dev_neg'];
%         SnewP = zeros(numel(metsListP),3);
%         SnewP(:,1) = coeffP;
%         SnewP(:,2) = coeffP;
%         SnewP(:,3) = -coeffP;
%         pointsModel = addMultipleReactions(pointsModel, {baseP,dposP,dnegP}, metsListP, SnewP, 'lb', [tgt,0,0], 'ub', [tgt,1000,10000]);

%         % newDietModel decomposition (for later analysis and consistency)
%         rxnIdxN = find(strcmp(newDietModel.rxns, roiName), 1);
%         if ~isempty(rxnIdxN)
%             stoichN = full(newDietModel.S(:, rxnIdxN));
%             metsMaskN = (stoichN ~= 0);
%             metsListN = newDietModel.mets(metsMaskN)';
%             coeffN = stoichN(metsMaskN);
%             evalc('[newDietModel,~,~] = removeRxns(newDietModel, {roiName});');
%             baseN = [roiName, '_base'];
%             dposN = [roiName, '_dev_pos'];
%             dnegN = [roiName, '_dev_neg'];
%             SnewN = zeros(numel(metsListN),3);
%             SnewN(:,1) = coeffN;
%             SnewN(:,2) = coeffN;
%             SnewN(:,3) = -coeffN;
            
%             newDietModel = addMultipleReactions(newDietModel, {baseN,dposN,dnegN}, metsListN, SnewN, 'lb', [tgt,0,0], 'ub', [tgt,1e6,1e6]);
%         end

%         % accumulate new rois and weights (duplicate weight for both dev reactions)
%         newRois{end+1} = dposP; 
%         newRois{end+1} = dnegP; 
%         if exist('roiWeights','var') && ~isempty(roiWeights) && length(roiWeights) >= i
%             newWeights = [newWeights, roiWeights(i), roiWeights(i)];
%         else
%             newWeights = [newWeights, 1, 1]; 
%         end
%     end
%     rois = newRois;
%     roiWeights = newWeights;
% end
disp("biomassFluxes")
disp(biomassFluxes)
disp("targetFlux")
disp(targetFlux)

if ~isempty(targetFlux)
    newRois = {};
    newWeights = [];
    for i = 1:length(rois)
        roiName = char(rois{i});
        tgt = targetFlux(i);        % control flux value
        pcos = biomassFluxes(i);    % PCOS flux value
        delta = pcos - tgt;         % deviation value

        % pointsModel decomposition
        rxnIdxP = find(strcmp(pointsModel.rxns, roiName), 1);
        if isempty(rxnIdxP)
            error(['ROI reaction not found in pointsModel: ', roiName]);
        end
        stoichP = full(pointsModel.S(:, rxnIdxP));
        metsMaskP = (stoichP ~= 0);
        metsListP = pointsModel.mets(metsMaskP)';
        coeffP = stoichP(metsMaskP);
        evalc('[pointsModel,~,~] = removeRxns(pointsModel, {roiName});');
        baseP = [roiName, '_base'];
        dposP = [roiName, '_dev_pos'];
        dnegP = [roiName, '_dev_neg'];
        SnewP = zeros(numel(metsListP), 3);
        SnewP(:, 1) = coeffP;           % base at control value
        SnewP(:, 2) = coeffP;           % dev_pos
        SnewP(:, 3) = -coeffP;          % dev_neg

        % Set bounds for deviation reactions
        if delta > 0
            devPosLB = delta;  devPosUB = delta;
            devNegLB = 0;      devNegUB = 0;
        elseif delta < 0
            devPosLB = 0;      devPosUB = 0;
            devNegLB = abs(delta); devNegUB = abs(delta);
        else
            devPosLB = 0;      devPosUB = 0;
            devNegLB = 0;      devNegUB = 0;
        end

        pointsModel = addMultipleReactions(pointsModel, ...
            {baseP, dposP, dnegP}, metsListP, SnewP, ...
            'lb', [tgt, devPosLB, devNegLB], ...
            'ub', [tgt, devPosUB, devNegUB]);

        % newDietModel decomposition (for later analysis and consistency)
        rxnIdxN = find(strcmp(newDietModel.rxns, roiName), 1);
        if ~isempty(rxnIdxN)
            stoichN = full(newDietModel.S(:, rxnIdxN));
            metsMaskN = (stoichN ~= 0);
            metsListN = newDietModel.mets(metsMaskN)';
            coeffN = stoichN(metsMaskN);
            evalc('[newDietModel,~,~] = removeRxns(newDietModel, {roiName});');
            baseN = [roiName, '_base'];
            dposN = [roiName, '_dev_pos'];
            dnegN = [roiName, '_dev_neg'];
            SnewN = zeros(numel(metsListN), 3);
            SnewN(:, 1) = coeffN;
            SnewN(:, 2) = coeffN;
            SnewN(:, 3) = -coeffN;

            if delta > 0
                devPosLB_N = delta;  devPosUB_N = delta;
                devNegLB_N = 0;      devNegUB_N = 0;
            elseif delta < 0
                devPosLB_N = 0;      devPosUB_N = 0;
                devNegLB_N = abs(delta); devNegUB_N = abs(delta);
            else
                devPosLB_N = 0;      devPosUB_N = 0;
                devNegLB_N = 0;      devNegUB_N = 0;
            end

            newDietModel = addMultipleReactions(newDietModel, ...
                {baseN, dposN, dnegN}, metsListN, SnewN, ...
                'lb', [tgt, devPosLB_N, devNegLB_N], ...
                'ub', [tgt, devPosUB_N, devNegUB_N]);
        end

        % accumulate new rois and weights (duplicate weight for both dev reactions)
        newRois{end + 1} = dposP;
        newRois{end + 1} = dnegP;
        if exist('roiWeights', 'var') && ~isempty(roiWeights) && length(roiWeights) >= i
            newWeights = [newWeights, roiWeights(i), roiWeights(i)];
        else
            newWeights = [newWeights, 1, 1];
        end
    end
    rois = newRois;
    roiWeights = newWeights;
end


% if ~isempty(objIndex) && any(model.ub(objIndex)~=model.lb(objIndex)) && all(~isnan(f1))
%     pointsModel.lb(objIndex)=OFS*f1;
% end

pointsModel=addMetabolite(pointsModel, 'unitOfFoodAdded[dP]');
pointsModel=addMetabolite(pointsModel, 'unitOfFoodRemoved[dP]');
pointsModel=addMetabolite(pointsModel, 'unitOfFoodChange[dP]');
pointsModel=addMetabolite(pointsModel, 'roiPoint[roiP]');
pointsModel=addMetabolite(pointsModel, 'point[P]');

%If necessary, add all diet exchange reactions to targetedDietRxns
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
disp(targetedDietRxns)



%Set up foodRemovalWeighting
if isempty(foodRemovalWeighting)
    %Add all Diet_EX reactions, and set weight to 1
    dietRxns=find(contains(pointsModel.rxns,'Diet_EX'));
    foodRemovalWeighting=[pointsModel.rxns(dietRxns),num2cell(ones(length(dietRxns),1))];
    fprintf('Food removal weighting not specified. All dietary reactions set to weight of one for removal.\n')
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




%Add "Food_Added" an "Food_Removed" reactions to points model
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





%Get roi Indexes
for i=1:length(rois)
    roiIndexP(i)=find(strcmp(pointsModel.rxns,rois{i}));
end
roiUB=pointsModel.ub(roiIndexP);
roiLB=pointsModel.lb(roiIndexP);

%replace roi function
stoich=pointsModel.S(:,roiIndexP);
if length(roiIndexP)>1
    metInd=find(any(stoich.'~=0));
    metsRoi=pointsModel.mets(any(stoich.'~=0)).';
    metsStoich=full(stoich(any(stoich.'~=0),:));
else
    metsRoi=pointsModel.mets(find(stoich~=0)).';
    metsStoich=full(stoich(find(stoich~=0)));
end

% weightVector = zeros(1,length(roiIndexP));
weightVector = ones(1, length(rois));
disp(weightVector);


metsCombined = [metsRoi, {'roiPoint[roiP]'}, {'point[P]'}]; % must be cell array for concatenation
nROI = length(rois);

% Initialize stoich matrix: rows = metabolites, cols = reactions (one reaction per ROI)
roiStoichMatrix = zeros(length(metsCombined), nROI);

% Assign weights to 'roiPoint[roiP]' metabolite (which is at position length(metsRoi)+1)
% and zero for 'point[P]' (last metabolite)
for i = 1:nROI
    roiStoichMatrix(length(metsRoi)+1, i) =  weightVector(i) * roiWeights(i);  
    % Negative because consumption in this reaction; check your convention 
    roiStoichMatrix(length(metsCombined), i) = 0;
end

% Add intermediate reactions: ROI -> roiPoint[roiP] (weighted)
pointsModel = addMultipleReactions(pointsModel, rois, metsCombined, roiStoichMatrix, ...
    'lb', roiLB, 'ub', roiUB);

% Add pooling reaction: roiPoint[roiP] -> point[P]
pointsModel = addMultipleReactions(pointsModel, {'Point_EX_roiPoints[roiP]_[P]'}, ...
    {'roiPoint[roiP]', 'point[P]'}, [-1; 1], 'lb', -1000000, 'ub', 1000000);




%Find solution
% pointsModel = changeObjective(pointsModel,'Point_EX_Point[P]');
% disp(pointsModel.rxns{pointsModel.c~=0});


% pointsModel.osenseStr = 'min';
% pointsModelSln = optimizeCbModel(pointsModel);

% --- Apply weighting INSIDE optimization (very important) ---
devIdx  = find(strcmp(pointsModel.rxns,'Point_EX_roiPoints[roiP]_[P]'));
foodIdx = find(strcmp(pointsModel.rxns,'Point_EX_unitOfFoodChange[dP]_[P]'));

Wdev  = 1;   % 80% weight for deviation
Wfood = 0;   % 20% weight for diet change

pointsModel.c(:) = 0;       % clear objective
pointsModel.c(devIdx)  = Wdev;
pointsModel.c(foodIdx) = Wfood;
pointsModel.osenseStr = 'min';

% --- NOW solve model (must be here) ---
pointsModelSln = optimizeCbModel(pointsModel);



% disp(pointsModelSln.v(find(strcmp(pointsModel.rxns,'Point_EX_roiPoints[roiP]_[P]'))));
% disp(pointsModelSln.v(find(strcmp(pointsModel.rxns,'Point_EX_Point[P]'))));
% disp('Checking fluxes through Food Added and Food Removed reactions:');

% for i = 1:length(foodAddedRxns)
%     rxnInd = find(strcmp(pointsModel.rxns, foodAddedRxns{i}), 1);
%     if ~isempty(rxnInd)
%         fprintf('%s flux = %.4f\n', foodAddedRxns{i}, pointsModelSln.v(rxnInd));
%     end
% end

% for i = 1:length(foodRemovedRxns)
%     rxnInd = find(strcmp(pointsModel.rxns, foodRemovedRxns{i}), 1);
%     if ~isempty(rxnInd)
%         fprintf('%s flux = %.4f\n', foodRemovedRxns{i}, pointsModelSln.v(rxnInd));
%     end
% end

% % Check pooling fluxes
% poolRxns = {'Point_EX_unitOfFoodAdded2Change[dp]', 'Point_EX_unitOfFoodRemoved2Change[dp]', 'Point_EX_unitOfFoodChange[dP]_[P]', 'Point_EX_Point[P]'};
% for i = 1:length(poolRxns)
%     rxnInd = find(strcmp(pointsModel.rxns, poolRxns{i}), 1);
%     if ~isempty(rxnInd)
%         fprintf('%s flux = %.4f\n', poolRxns{i}, pointsModelSln.v(rxnInd));
%     else
%         fprintf('Pooling reaction %s not found in model\n', poolRxns{i});
%     end
% end

% for i = 1:length(foodAddedRxns)
%     rxnID = find(strcmp(pointsModel.rxns, foodAddedRxns{i}));
%     if ~isempty(rxnID)
%         fprintf('%s bounds: LB=%.2f, UB=%.2f\n', foodAddedRxns{i}, pointsModel.lb(rxnID), pointsModel.ub(rxnID));
%     end
% end



biomassRxns = {'SK_biomass_maintenance', 'AD_biomass_maintenance', ...
               'GN_biomass_maintenance', 'OO_biomass_maintenance', ...
               'EN_biomass_maintenance',  'SK_PGESr_base', ...
               'GN_PGSr_base', 'GN_HMR_2581_base', 'OO_HMR_0987_base', ...
               'OO_LTC4CP_base', 'GN_RE1796R_base', 'EN_HAS2_base', ...
               'AD_ARGSL_base', 'AD_SPHK11_base', 'SK_PGESr_dev_pos', ...
               'SK_PGESr_dev_neg', 'GN_PGSr_dev_pos', 'GN_PGSr_dev_neg', ...
               'GN_HMR_2581_dev_pos', 'GN_HMR_2581_dev_neg', 'OO_HMR_0987_dev_pos', ...
               'OO_HMR_0987_dev_neg', 'OO_LTC4CP_dev_pos', 'OO_LTC4CP_dev_neg', ...
               'GN_RE1796R_dev_pos', 'GN_RE1796R_dev_neg', 'EN_HAS2_dev_pos', ...
               'EN_HAS2_dev_neg', 'AD_ARGSL_dev_pos', 'AD_ARGSL_dev_neg', ...
               'AD_SPHK11_dev_pos', 'AD_SPHK11_dev_neg'};

% [isPresent, rxnIndices] = ismember(biomassRxns, pointsModel.rxns);

% fprintf('\nFlux values for ROI reactions after minimization:\n');
% for i = 1:length(biomassRxns)
%     if isPresent(i)
%         fprintf('%s = %.4f\n', biomassRxns{i}, pointsModelSln.v(rxnIndices(i)));
%         fprintf('Upper bound for %s = %.4f\n', biomassRxns{i}, pointsModel.ub(rxnIndices(i)));
%     else
%         fprintf('%s not found in pointsModel\n', biomassRxns{i});
%     end
% end

if strcmp(display,'on')
    disp(['Solution points =',num2str(pointsModelSln.f)])
end

if strcmp(display, 'on')
    disp([num2str(pointsModelSln.v(find(strcmp(pointsModel.rxns,'Point_EX_unitOfFoodChange[dP]_[P]')))),' come from diet']);
    disp([num2str(pointsModelSln.v(find(strcmp(pointsModel.rxns,'Point_EX_roiPoints[roiP]_[P]')))),' come from rois']);
end
% --- Weighting system: 80% deviation, 20% food ---
dietFlux = pointsModelSln.v(find(strcmp(pointsModel.rxns,'Point_EX_unitOfFoodChange[dP]_[P]')));
devFlux  = pointsModelSln.v(find(strcmp(pointsModel.rxns,'Point_EX_roiPoints[roiP]_[P]')));

Wdev  = 0.80; % 80% weight for deviation
Wfood = 0.20; % 20% weight for food

weightedTotalPoints = (Wdev * devFlux) + (Wfood * dietFlux);

if strcmp(display,'on')
    disp(['Weighted Total Points = ', num2str(weightedTotalPoints)]);
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
devCount = 0;
for i=1:length(rois)
    ind=find(strcmp(pointsModel.rxns,rois{i}));
    fluxVal = pointsModelSln.v(ind(1));
    isDev = contains(rois{i},'_dev_pos') || contains(rois{i},'_dev_neg');
    if isDev && fluxVal == 0
        continue; % skip printing zero deviation fluxes
    end
    if isDev && fluxVal ~= 0
        devCount = devCount + 1;
    end
    if strcmp(display,'on')
        disp(['   ',rois{i},' flux = ', num2str(fluxVal)])
    end
    roiFlux(i)=fluxVal;
end
if strcmp(display,'on')
    disp(['Nonzero deviation flux reactions: ', num2str(devCount)])
end
if strcmp(slnType,'Quick')
    detailedAnalysis=[];
    return
end


%Find new obj flux with new diet

if isempty(targetFlux)
    if strcmp(display,'on')
        disp('Detailed Analysis:')
    end
    if ~isempty(objIndex) && any(model.ub(objIndex)~=model.lb(objIndex))
        model_Obj = optimizeCbModel(newDietModel);
        f2 = model_Obj.f;
        % Constrain all tissue biomass maintenance reactions
        tissueBiomassRxns = {'SK_biomass_maintenance','AD_biomass_maintenance','GN_biomass_maintenance','OO_biomass_maintenance','EN_biomass_maintenance'};
        weights = [0.2, 0.2, 0.2, 0.2, 0.2];
        for i = 1:length(tissueBiomassRxns)
            rxnIdx = find(strcmp(newDietModel.rxns, tissueBiomassRxns{i}));
            if ~isempty(rxnIdx)
                newDietModel.c(rxnIdx) = weights(i);
                newDietModel = changeRxnBounds(newDietModel, tissueBiomassRxns{i}, model.ub(rxnIdx), 'u');
                newDietModel = changeRxnBounds(newDietModel, tissueBiomassRxns{i}, model.lb(rxnIdx), 'l');
            end
        end
        if strcmp(display,'on')
            disp([' ','Original objective max flux =',num2str(f1), ' & New objective max flux =', num2str(f2)])
        end
    else
        f2 = model.ub(objIndex);
        if strcmp(display,'on')
            disp([' ','Original objective max flux =',num2str(f1), ' & New objective max flux =', num2str(f2)])
        end
    end

    % Compute new min max ranges for roi with new diet
    for i=1:length(rois)
        if strcmp(display,'on')
            disp([' ' rois{i}])
        end
        tmpModel = newDietModel;
        tmpModel.c(:) = 0;
        rxnIdx = find(strcmp(tmpModel.rxns, rois{i}));
        if ~isempty(rxnIdx)
            tmpModel.c(rxnIdx) = 1;
        end
        tmpModel.osenseStr = 'min';
        tmp = optimizeCbModel(tmpModel);
        detailedAnalysis.(['Rxn',num2str(i)]).min.ND = tmp;
        NroiFluxMin(i)=tmp.v(roiIndexO(i));
        tmpModel.osenseStr = 'max';
        tmp = optimizeCbModel(tmpModel);
        detailedAnalysis.(['Rxn',num2str(i)]).max.ND = tmp;
        NroiFluxMax(i)=tmp.v(roiIndexO(i));
        if strcmp(display,'on')
            disp(['   Original Diet RoI range = ', num2str(OroiFluxMin(i)), ':', num2str(OroiFluxMax(i))])
            disp(['   New Diet RoI range = ', num2str(NroiFluxMin(i)), ':', num2str(NroiFluxMax(i))])
        end
    end
end

newDietModel.ub(objIndex)=model.ub(objIndex);
newDietModel.lb(objIndex)=model.lb(objIndex);
% Set objective as weighted sum over all tissue biomass maintenance reactions
tissueBiomassRxns = {'SK_biomass_maintenance','AD_biomass_maintenance','GN_biomass_maintenance','OO_biomass_maintenance','EN_biomass_maintenance'};
weights = [0.2, 0.2, 0.2, 0.2, 0.2];
newDietModel.c(:) = 0;
for i = 1:length(tissueBiomassRxns)
    rxnIdx = find(strcmp(newDietModel.rxns, tissueBiomassRxns{i}));
    if ~isempty(rxnIdx)
        newDietModel.c(rxnIdx) = weights(i);
    end
end
newDietModel.osenseStr = objMinMax;

if strcmp(display,'on')
    disp('_____________________________________________________')
end

% ==========================================
% CLEAN newDietModel: KEEP base, REMOVE dev, RESTORE names
% ==========================================

rxnsToRemove = {};

% Remove ONLY deviation reactions (NOT base)
rxnsToRemove = [rxnsToRemove; newDietModel.rxns(contains(newDietModel.rxns, '_dev_pos'))];
rxnsToRemove = [rxnsToRemove; newDietModel.rxns(contains(newDietModel.rxns, '_dev_neg'))];

% Remove points & food optimization reactions
rxnsToRemove = [rxnsToRemove; newDietModel.rxns(contains(newDietModel.rxns, 'Point_EX_'))];
rxnsToRemove = [rxnsToRemove; newDietModel.rxns(contains(newDietModel.rxns, 'Food_Added_EX_'))];
rxnsToRemove = [rxnsToRemove; newDietModel.rxns(contains(newDietModel.rxns, 'Food_Removed_EX_'))];

rxnsToRemove = unique(rxnsToRemove);

% Remove helper reactions
if ~isempty(rxnsToRemove)
    evalc('[newDietModel,~,~] = removeRxns(newDietModel, rxnsToRemove);');
end

% Rename _base reactions back to original names
baseRxns = newDietModel.rxns(contains(newDietModel.rxns, '_base'));

for i = 1:length(baseRxns)
    oldName = baseRxns{i};
    newName = erase(oldName, '_base');
    newDietModel.rxns(strcmp(newDietModel.rxns, oldName)) = {newName};
end

disp(['Removed ', num2str(length(rxnsToRemove)), ' helper reactions']);
disp(['Restored ', num2str(length(baseRxns)), ' base reactions to original names']);


end
