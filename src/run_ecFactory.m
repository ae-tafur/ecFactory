function [optStrain,remaining] = run_ecFactory(model,targetRxn,csRxn,csMW,expYield,essential,outputPath,graphPlot,modelAdapter)
% run_ecFactory
%   Run ecFactory methods to fin and optimal mutant based on ecFSEOF
%   results for a specfified production target.
%
% Input:
%   ecModel         an ecModel in GECKO 3 format (with ecModel.ec structure).
%   targetRxn       rxn ID for the production target reaction, a exchange
%                   reaction is recommended.
%   csRxn           rxn ID for the main carbon source uptake reaction.
%   csMW            Molecular weigth for the carbon source in grams/mmol.
%   nSteps          a suboptimal yield value to run ecFactory (e.g., 0.1 g/g)
%   essential       file name for essential list of genes. (Optional,
%                   default ecFactory will look the 'data' sub-folder from
%                   GECKO v3 structure and taken from
%                   modelAdapter, e.g. data/essential_genes.tsv)
%   outputPath      path to save results output files.(Optional, default
%                   in the 'output' sub-folder taken from
%                   modelAdapter, e.g. output/ecFSEOF_genes.tsv
%   file_rxns       true or false if a plot with targets interaction at
%                   metabolite level is desired. (Optional, default false)
%   modelAdapter    a loaded model adapter. (Optional, will otherwise use
%                   the default model adapter)
%
% Output:
%   optStrain       Optimal strain
%   remaining       gene target suggested for optimal strain
%
% Usage:
%   [optStrain,remaining] = run_ecFactory(model,targetRxn,csRxn,csMW,expYield)

if nargin < 9 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 8 || isempty(graphPlot)
    graphPlot = false;
end

if nargin < 7 || isempty(outputPath)
    outputPath = fullfile(params.path,'output');
end

if nargin < 6 || isempty(essential)
    % Read file with essential genes list
    essential = readtable(fullfile(params.path,'data','essential_genes.txt'),'Delimiter','\t');
end

if ~exist(outputPath, 'dir')
    mkdir(outputPath)
end

essential = strtrim(essential.Ids);

% Parameters for FSEOF method
nSteps     = 16; % number of FBA steps in ecFSEOF
alphaLims  = [0.5*expYield 2*expYield]; % biomass yield limits for ecFSEOF
file_genes = fullfile(outputPath,'ecFSEOF_genes.tsv');
file_rxns  = fullfile(outputPath,'ecFSEOF_rxns.tsv');

% Method parameters
tol  = 1E-8;  %numeric tolerance for determining non-zero enzyme usages
OEF  = 2;     %overexpression factor for enzyme targets
KDF  = 0.5;   %down-regulation factor for enzyme targets
thresholds = [0.5 1.05]; %K-score thresholds for valid gene targets
delLimit   = 0.05; %K-score limit for considering a target as deletion

% Get relevant rxn indexes
targetRxnIdx = getIndexes(model, targetRxn,'rxns');
csRxnIdx = getIndexes(model, csRxn,'rxns');
bioRxnIdx = getIndexes(model, params.bioRxn,'rxns');

% perform FBA to determine uptake rate
sol = solveLP(model);
uptake = abs(sol.x(csRxnIdx));

% ecModel verification steps
% model = check_enzyme_fields(model); %% add mean MW to those that have
% this missing
if ~isempty(targetRxnIdx)
    %Check if model can carry flux for the target rxn
    [~, maxFlux, ~] = getAllowedBounds(model, targetRxn);
    if maxFlux
        disp(['* Your ecModel can carry flux through the reaction: ' model.rxnNames{targetRxnIdx}])
    else
        disp(['* Your ecModel cannot carry flux through the reaction: ' model.rxnNames{targetRxnIdx} ', please check the applied constraints'])
    end
else
    error('The provided target reaction is not part of the ecModel.rxns field')
end


% 1.- Run FSEOF to find gene candidates

fprintf('\n1.- **** Running ecFSEOF method (from GECKO utilities) **** \n\n')
results = run_ecFSEOF(model,targetRxn,csRxn,alphaLims,nSteps,file_genes,file_rxns);
genes   = results.geneTable(:,1);
fprintf(['\n  ecFSEOF returned ' num2str(length(genes)) ' targets \n'])

%Format results table
geneShorts = results.geneTable(:,2);
k_scores   = cell2mat(results.geneTable(:,3));
actions    = cell(numel(k_scores),1);
actions(k_scores >= thresholds(2)) = {'OE'};
actions(k_scores < thresholds(1))  = {'KD'};
actions(k_scores < delLimit) = {'KO'};

%Identify candidate genes in model enzymes
fprintf('\n  Extracting enzymatic information for target genes \n')
[~,iB] = ismember(genes,model.ec.genes);
MWeigths = nan(numel(iB),1);
candidates = cell(numel(iB),1);
pathways   = cell(numel(iB),1);
%optimize!
for i=1:numel(iB)
    if iB(i)>0
        candidates(i) = model.ec.enzymes(iB(i));
        MWeigths(i)   = model.ec.mw(iB(i));
        % pathways(i)   = model.subSystems(iB(i)); %% revisar
    end
end

%Get results table
candidates = table(genes,candidates,geneShorts,MWeigths,pathways,actions,k_scores, ...
    'VariableNames',{'genes' 'enzymes' 'shortNames' 'MWs' 'pathways' 'actions' 'k_scores'});
%Keep results that comply with the specified K-score thresholds
fprintf(['\n  Removing targets ' num2str(thresholds(1)) ' < K_score < ' num2str(thresholds(2)) '\n'])
toKeep     = find((candidates.k_scores>=thresholds(2)|candidates.k_scores<=thresholds(1)));
candidates = candidates(toKeep,:);
fprintf(['\n    * ' num2str(height(candidates)) ' gene targets remain \n']) 

%2.- Add flux leak targets (those genes not optimal for production that may
%consume the product of interest. (probaly extend the approach to inmediate
%precurssors)
fprintf('\n2.- **** Find flux leak targets to block **** \n\n')
candidates = find_flux_leaks(candidates,targetRxnIdx,model);
fprintf(['\n    * ' num2str(height(candidates)) ' gene targets remain \n']) 

% 3.- discard essential genes from deletion targets
fprintf('\n3.- **** Removing essential genes from KD and KO targets list **** \n')
[~,iB]    = ismember(candidates.genes,essential);
toRemove  = iB & candidates.k_scores<=delLimit;
candidates(toRemove,:) = [];
fprintf(['\n    * ' num2str(height(candidates)) ' gene targets remain \n'])
writetable(candidates,fullfile(outputPath, 'candidates_L1.txt'),'Delimiter','\t','QuoteStrings',false);
proteins = strcat('usage_prot_',candidates.enzymes);
[~,enz_pos] = ismember(proteins,model.rxns);
candidates.enz_pos = enz_pos;

% 4.- Construct Genes-metabolites network for classification of targets
fprintf('\n4.- **** Construct Genes-metabolites network for classification of targets **** \n')
%Get Genes-metabolites network
fprintf('\n  Constructing genes-metabolites graph \n')
[GeneMetMatrix,~,Gconect] = getMetGeneMatrix(model,candidates.genes);
%Get independent genes from GeneMetMatrix
fprintf('\n  Obtain redundant vectors in genes-metabolites graph (redundant targets) \n')

%generate gene-gene matrix and identify linearly independent targets
[indGenes,G2Gmatrix,~] = getGeneDepMatrix(GeneMetMatrix);
%find unique targets (with no isoenzymes or not part of complexes)
candidates.unique = indGenes;
%number of metabolites connected to each gene
candidates.conectivity = Gconect.mets_number; 
%Get gene target groups (those connected to exactly the same metabolites)
[~,groups]        = getGenesGroups(G2Gmatrix,candidates.genes);
candidates.groups = groups;

% 5.- Enzyme usage variability analysis (EVA) and prioritization of targets
fprintf('\n5.- **** Running EUVA for optimal production conditions (minimal biomass) **** \n')
tempModel = model;
fprintf(['\n    * Fixed ' erase(model.rxnNames{csRxnIdx}, ' exchange') ' uptake rate to ' num2str(uptake) '\n'])
%Fix unit C source uptake
tempModel = setParam(tempModel, 'var', csRxn, -uptake, tol);
fprintf('    * Fixed suboptimal biomass production, according to provided experimental yield.')
%Fix suboptimal experimental biomass yield conditions
V_bio = expYield*csMW;
tempModel.lb(bioRxnIdx) = V_bio;
fprintf([' V_bio = ' num2str(V_bio) ' h-1 \n'])
fprintf('    * Production rate constrained to its maximum attainable value.')
%Get and fix optimal production rate
tempModel = setParam(tempModel, 'obj', targetRxn, +1);
sol       = solveLP(tempModel,1);
max_prod  = sol.x(targetRxnIdx);
tempModel = setParam(tempModel, 'var', targetRxn, max_prod, tol);
fprintf([' V_prod_max = ' num2str(max_prod) ' mmol/gDW/h \n'])

%Run FVA for all enzyme usages subject to fixed CUR and Grates
EVAtable = enzymeUsage_FVA(tempModel,candidates.enzymes);

fprintf('\n  Discard targets according to EUVA results for optimal production **** \n')
%Classify targets according to enzyme variability ranges
candidateUsages = EVAtable.enzUsage;
candidates.EV_type = cell(height(candidates),1);
candidates.EV_type(:) = {''};
idxs = find(EVAtable.minUsage<=tol & EVAtable.maxUsage<=tol);
candidates.EV_type(idxs) = {'unusable'};
idxs = find(EVAtable.minUsage>tol);
candidates.EV_type(idxs) = {'essential'}; %enzymes needed for optimal production
idxs = find(EVAtable.minUsage<=tol & EVAtable.maxUsage>tol);
candidates.EV_type(idxs) = {'totally_variable'}; %enzymes that can take any flux (ussually isoenzymes)
idxs = find((EVAtable.minUsage./EVAtable.maxUsage)>=0.99 & EVAtable.minUsage>tol);
candidates.EV_type(idxs) = {'essential_tightly_const'}; %enzymes that are needed up to its maximum capacity
idxs = find(strcmpi(candidates.EV_type,'totally_variable') & candidateUsages>tol);
candidates.EV_type(idxs) = {'production_opt'};  %isoenzymes that are chosen for optimal production 
idxs = find(strcmpi(candidates.EV_type,'production_opt') & (candidateUsages./EVAtable.maxUsage)>=0.99);
candidates.EV_type(idxs) = {'production_opt_tight'}; %isoenzymes that are chosen for optimal production, up to its maximum capacity
idxs = find(strcmpi(candidates.EV_type,'totally_variable') & candidateUsages<=tol);
candidates.EV_type(idxs) = {'unnecessary_prod'}; %isoenzymes that are not chosen for optimal production
%Append EUVA results to target candidates table
candidates.OE(strcmpi(candidates.actions,'OE')) = OEF;
candidates.OE(strcmpi(candidates.actions,'KD')) = KDF;
candidates.OE(strcmpi(candidates.actions,'KO')) = 0;
candidates.minUsage = EVAtable.minUsage;
candidates.maxUsage = EVAtable.maxUsage;
candidates.enzUsage = candidateUsages;
%Discard enzymes 
fprintf('    * Discard OE targets with lb=ub=0 \n')
toRemove = (strcmpi(candidates.EV_type,'unusable') & strcmpi(candidates.actions,'OE')) | isnan(candidates.minUsage);
candidates(toRemove,:) = [];
fprintf('    * Discard enzymes essential for production from deletion targets \n')
toRemove = (strcmpi(candidates.EV_type,'essential_tightly_const') | strcmpi(candidates.EV_type,'essential')) & ...
       (candidates.k_scores<=delLimit);
candidates(toRemove,:) = [];       
fprintf('    * Discard isoenzyme groups for KD/KO that contain an optimal isoform (redundant groups that increase mutant complexity) \n')
toRemove = [];
for k=1:max(candidates.groups)
    idx = find(candidates.groups==k & ~strcmpi(candidates.actions,'OE'));
    if length(idx)>1
        if ~isempty(candidates.enzymes(idx(1))) & any(candidates.enzUsage(idx)>0)
            toRemove = [toRemove;idx];
        end
    end
end
candidates(toRemove,:) = [];
fprintf(['\n    * ' num2str(height(candidates)) ' gene targets remain \n'])

%6.- EUVA for suboptimal biomasss production subject to a minimal (1%)
% production rate of the target product and a unit CS uptake rate
%Get max biomass 
fprintf('\n6.-  **** Running EUVA for reference strain  **** \n')
fprintf(['    * Fixed ' erase(model.rxnNames{csRxnIdx}, ' exchange') ' uptake rate to ' num2str(uptake) '\n'])
fprintf('    * Production rate subject to a LB of 1%% of its max. value \n')
fprintf('    * Biomass production fixed to its maximum attainable value \n')
tempModel = setParam(tempModel, 'obj', params.bioRxn, +1);
tempModel.lb(targetRxnIdx) = 0.01*max_prod;
tempModel.ub(targetRxnIdx) = 1000;
sol       = solveLP(tempModel,1);
maxVBio   = sol.x(bioRxnIdx);
%Fix optimal biomass formation rate
tempModel = setParam(tempModel, 'var', params.bioRxn, maxVBio, tol);
%run EUVA for optimal biomass formation
EVAbio = enzymeUsage_FVA(tempModel,candidates.enzymes);
candidates.minUsageBio = EVAbio.minUsage;
candidates.maxUsageBio = EVAbio.maxUsage;
candidates.enzUsageBio = EVAbio.enzUsage;

%discard something 
fprintf('\n  Discard targets according to EUVA results for reference strain **** \n')
candidates = compare_EUVR(candidates);
fprintf('    * Discarding enzymes with inconsistent enzyme usage variability patterns')
toRemove = strcmpi(candidates.EUV_comparison,'embedded') | ...
           (~strcmpi(candidates.actions,'OE') & contains(candidates.EUV_comparison,'up_')) | ...
           (strcmpi(candidates.actions,'OE') & contains(candidates.EUV_comparison,'down_'));% |...
%disp(candidates(toRemove,:))
candidates(toRemove,:) = [];
fprintf(['\n    * ' num2str(height(candidates)) ' gene targets remain \n'])

%7.-
fprintf('\n7.-  **** Rank targets by priority levels **** \n\n')
%Rank targets by priority
candidates.priority = zeros(height(candidates),1)+3;
candidates.priority(contains(candidates.EUV_comparison,'up_') | contains(candidates.EUV_comparison,'down_')) = 2; % second priority for genes with overlaped demand ranges
candidates.priority(contains(candidates.EUV_comparison,'_distinct')) = 1; %higher priority to the clearly up-down regulated genes
fprintf(['    * ' num2str(sum(candidates.priority==1)) ' targets with priority level 1 (distinct enzyme demand levels between optimal production and optimal biomass cases) \n']) 
fprintf(['    * ' num2str(sum(candidates.priority==2)) ' targets with priority level 2 (overlapped enzyme demand levels between optimal production and optimal biomass cases) \n']) 
fprintf(['    * ' num2str(sum(candidates.priority==3)) ' targets with priority level 3 (other cases) \n']) 

%Rename enzyme variability type for KDs and KOs according to their
%variability ranges for maximum biomass production
fprintf('\n Classifying targets according to their enzyme usage variability range type \n')
idxs = ~strcmpi(candidates.actions,'OE') & candidates.enzUsage < candidates.enzUsageBio & candidates.enzUsageBio>0 & ~contains(candidates.EV_type,'unnecessary');
candidates.EV_type(idxs) = {'biomass_opt'};
bioRatio = V_bio/maxVBio;
ratios   = candidates.enzUsage./candidates.enzUsageBio;
idxs     = ratios < bioRatio+1E-9 & ratios > bioRatio-1E-9;
candidates.EV_type(idxs) = {'biomass_coupled'};
writetable(candidates,fullfile(outputPath, 'candidates_L2.txt'),'Delimiter','\t','QuoteStrings',false);

% 8.- Combine targets
fprintf('\n8.-  **** Find an optimal combination of remaining targets **** \n')
%Unconstrain CUR, unconstrain product formation 
%and set a minimal biomass formation
tempModel = setParam(tempModel,'ub',csRxn,0);
tempModel = setParam(tempModel,'lb',csRxn,-1000);
tempModel = setParam(tempModel,'ub',params.bioRxn,1000);
tempModel = setParam(tempModel,'lb',params.bioRxn,0.99*V_bio);
tempModel = setParam(tempModel,'ub',targetRxn,1000);
tempModel = setParam(tempModel,'lb',targetRxn,0);
%set Max product formation as objective function
tempModel = setParam(tempModel,'obj',targetRxn,+1);
%constrain enzyme usages with optimal biomass formation profile (WT)
toChange = isnan(candidates.maxUsageBio);
tempModel.lb(candidates.enz_pos(~toChange)) = 1.01*-candidates.maxUsageBio(~toChange);
tempModel.lb(candidates.enz_pos(toChange))  = 1.01*-candidates.enzUsageBio(toChange);
tempModel.ub(candidates.enz_pos)            = 0.99*-candidates.enzUsageBio;
%Run mechanistic validation of targets
[optStrain,remaining] = constructMinimalMutant(tempModel,candidates,targetRxn,csRxn,csMW,params.bioRxn);
%get a gene-mets graph with the remaining targets 
if graphPlot
    MetsIndxs    = (~contains(optStrain.metNames,'prot_'));
    nodeMets     = optStrain.mets(MetsIndxs);
    toKeep       = (remaining.priority>0 & remaining.conectivity<15);
    [GeneMetMatrix,~,~] = getMetGeneMatrix(optStrain,remaining.genes);
    tempGMmatrix = GeneMetMatrix(:,toKeep);
    getMGgraph(tempGMmatrix,nodeMets,tempModel,'force',remaining(toKeep,:));
end
fprintf(['\n    * The predicted optimal strain contains ' num2str(height(remaining)) ' gene modifications \n']) 

writetable(remaining,fullfile(outputPath, 'candidates_L3.txt'),'Delimiter','\t','QuoteStrings',false);

% 9. Generate transporter targets file (lists a number of transport steps
% with no enzymatic annotation that are relevant for enhancing target
% product formation.

fprintf('\n9.-  **** Find transporter reactions with no enzyme association predicted as targets by ecFSEOF **** \n')
transpTargets = getTransportTargets(cell2table(results.rxnsTable, ...
    'VariableNames', {'rxn_IDs' 'rxnNames' 'K_score' 'grRules' 'rxn_formula'}),tempModel);
fprintf(['\n    * ' num2str(height(transpTargets)) ' transport reaction targets were found \n']) 
writetable(transpTargets,fullfile(outputPath, 'transporter_targets.txt'),'Delimiter','\t','QuoteStrings',false);
end