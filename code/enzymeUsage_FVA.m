function FVAtable = enzymeUsage_FVA(model,enzymes)
% enzymeUsage_FVA
%   Perf
%
%   Usage: FVAtable = enzymeUsage_FVA(model,enzymes)
%
%   Ivan Domenzain.     Last edited 2020-05-27

if nargin < 2
    enzymes = model.ec.enzymes;
end

%Get parsimonious protein usages
model = setParam(model, 'obj', 'prot_pool_exchange', 1); % Minimize protein pool
solUsgs = solveLP(model);
%poolIdx = find(strcmpi(model.rxns,'prot_pool_exchange'));
%model = setParam(model, 'lb', 'prot_pool_exchange', solUsgs.x(poolIdx));

%initialize variables
enzIDs    = cell(length(enzymes),1);
enzpUsgs  = zeros(length(enzymes),1);
minUsgs   = zeros(length(enzymes),1);
maxUsgs   = zeros(length(enzymes),1);
ranges    = zeros(length(enzymes),1);

if ~isempty(solUsgs.x)
    % pUsgs = sol.x(prot_indxs);
    % Loop through all the provided enzymes
    for i=1:length(enzymes)
        if ~isempty(enzymes{i})
            rxnIndx = find(strcmpi(model.rxns, strcat('usage_prot_', enzymes{i})));
            enzpUsgs(i) = solUsgs.x(rxnIndx);
            % enzIndx = find(strcmpi(model.ec.enzymes,enzymes{i}),1);
            enzIDs(i)  = enzymes(i);
            if ~isempty(rxnIndx)
                model  = setParam(model, 'obj', rxnIndx, +1);
                sol    = solveLP(model);
                if ~isempty(sol.f)
                    minUsgs(i) = sol.x(rxnIndx);
                    model   = setParam(model, 'obj', rxnIndx, -1);
                    sol     = solveLP(model);
                    if ~isempty(sol.f)
                        %disp(['Ready with enzyme #' num2str(i) ' ' model.enzymes{enzIndx}])
                        maxUsgs(i) = sol.x(rxnIndx);
                    end
                end
                ranges(i) = (minUsgs(i)-maxUsgs(i));
            end
        else
            ranges(i)    = NaN;
            minUsgs(i)   = NaN;
            maxUsgs(i)   = NaN;
            enzpUsgs(i)  = NaN;
            %enzIDs{i}    = 'empty';
        end
    end
else
    disp('Model is not feasible')
end
varNamesT = {'enzymes' 'ranges' 'minU' 'maxU' 'pU'};
FVAtable  = table(enzIDs,ranges,minUsgs,maxUsgs,enzpUsgs,'VariableNames', varNamesT);
end