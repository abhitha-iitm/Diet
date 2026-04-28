baseFolder = 'Diet_models';
controlFolder = fullfile(baseFolder,'control');
pcosFolder = fullfile(baseFolder,'pcos');

if ~exist(baseFolder,'dir')
    mkdir(baseFolder);
end

if ~exist(controlFolder,'dir')
    mkdir(controlFolder);
end

if ~exist(pcosFolder,'dir')
    mkdir(pcosFolder);
end



model_CT = readCbModel('./Multitissue_Models/MulModel_CT.mat');
model_PC = readCbModel('./Multitissue_Models/MulModel_PC.mat');

model_CT = convert_EX_to_diet(model_CT);
model_PC = convert_EX_to_diet(model_PC);


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
};


for d = 1:length(dietFiles)

    fprintf('Processing diet: %s\n', dietFiles{d});

    dietPath = fullfile('Diets',dietFiles{d});
    dietName = erase(dietFiles{d},'.tsv');

    %% Apply diet to CONTROL model
    [model_CT_diet,~,~] = setDietBoundsFromFile(model_CT,dietPath);

    %% Apply diet to PCOS model
    [model_PC_diet,~,~] = setDietBoundsFromFile(model_PC,dietPath);


    %% Save CONTROL model
    save(fullfile(controlFolder,...
        ['CT_' dietName '.mat']),...
        'model_CT_diet');


    %% Save PCOS model
    save(fullfile(pcosFolder,...
        ['PC_' dietName '.mat']),...
        'model_PC_diet');

end