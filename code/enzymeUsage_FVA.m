function fvaTable = enzymeUsage_FVA(model,enzymes)
% enzymeUsage_FVA
%   Perf
%
%   Usage: FVAtable = enzymeUsage_FVA(model,enzymes)
%
%   Ivan Domenzain.     Last edited 2020-05-27

if nargin < 2
    enzymes = model.ec.enzymes;
end

% Get protein usage reaction index 
[~,rxnIdx] = ismember(strcat('usage_prot_', enzymes), model.rxns);

% Initialize variables
enzIds  = erase(model.rxns(rxnIdx),'usage_prot_'); % Only keep members

% Get parsimonious protein usages
model = setParam(model, 'obj', 'prot_pool_exchange', 1); % Minimize protein pool
sol = solveLP(model, 1);

if ~isempty(sol.x)

    % Protein usage
    enzUsage = sol.x(rxnIdx);

    % Get maximum and minimum protein usage
    [maxUsage, minUsage, ~] = getAllowedBounds(model, rxnIdx);

    % Use positve numbers
    enzUsage = abs(enzUsage);
    minUsage = abs(minUsage);
    maxUsage = abs(maxUsage);

    % Set values which are under solver tolerance
    enzUsage(enzUsage < 1e-8) = 0;
    minUsage(minUsage < 1e-8) = 0;
    maxUsage(maxUsage < 1e-8) = 1e-8;

    ranges = maxUsage - minUsage;

else
    disp('Model is not feasible')
end
varNamesT = {'enzymes' 'ranges' 'minUsage' 'maxUsage' 'enzUsage'};
fvaTable  = table(enzIds,ranges,minUsage,maxUsage,enzUsage,'VariableNames', varNamesT);
end