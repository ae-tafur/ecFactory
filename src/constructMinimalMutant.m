function [optMutant,remaining] = constructMinimalMutant(model,candidates,targetRxn,csRxn,csMW,bioRxn)
tol =-1E-12;%-1E-6;
tempModel = model;

targetRxnIdx = getIndexes(model, targetRxn,'rxns');
csRxnIdx = getIndexes(model, csRxn,'rxns');
bioRxnIdx = getIndexes(model, bioRxn,'rxns');

%Get max WT production rate
WTsol_prod = solveECmodel(tempModel,tempModel,'pFBA','prot_pool_exchange');
%Get max WT production yield
WTsol_yield = solveECmodel(tempModel,tempModel,'pFBA',csRxn);
WTprodR = WTsol_prod(targetRxnIdx);
WTyield = WTsol_yield(targetRxnIdx)/abs(WTsol_yield(csRxnIdx));
%get mutant with all modifications
optMutant=model;
toRemove = [];
[candidates,~] = sortrows(candidates,{'priority' 'k_scores'},{'ascend' 'ascend'});

for i=1:height(candidates)
    gene   = candidates.genes(i);
    action = candidates.actions(i);
    mutF   = 1;%candidates.OE(i);
    enzIdx = candidates.enz_pos(i);
    tempModel = optMutant;

    if enzIdx>0
        % Validate if it was posible to get maximum values
        if ~isnan(candidates.maxUsage(i))
            tempModel.lb(enzIdx) = 1.01*-candidates.maxUsage(i);
        % Otherwise use the enzyme usage value
        else
            tempModel.lb(enzIdx) = 1.01*-candidates.enzUsage(i);
        end
        tempModel.ub(enzIdx) = 0;
        if tempModel.ub(enzIdx) <=tempModel.lb(enzIdx)
            % tempModel.lb(enzIdx) = 0.95*tempModel.ub(enzIdx);
        end

        optMutant = tempModel;
        %for reactions without enzymatic reaction
    else
        if strcmpi(action,'OE')
            enzUsage = 1.01*candidates.maxUsage(i);
            if enzUsage <= 1E-15
                enzUsage = 1.01*candidates.maxUsage(i);
            end
        else
            enzUsage = candidates.pUsage(i);
            if strcmpi(action,'KO')
                action = {'KD'};
                enzUsage = 0;
            end
        end
        modifications = {gene action mutF};

        [tempModel,success] = getMutantModel(optMutant,modifications,enzUsage);
        if success
            optMutant = tempModel;
        else
            toRemove = [toRemove; i];
        end
    end
end

if ~isempty(toRemove)
    disp('The following gene modifications are not compatible with the rest of remaining candidate targets')
    disp(candidates(toRemove,[1 2 3 6]))
    candidates(toRemove,:) = [];
end

optMutant = setParam(optMutant,'obj',targetRxn,1);
%obtain optimal production rate and yield
[mutSol_r,~] = solveECmodel(optMutant,model,'pFBA','prot_pool_exchange');
[mutSol_y,~] = solveECmodel(optMutant,model,'pFBA',csRxn);
OptprodR = mutSol_y(targetRxnIdx);
Optyield = mutSol_y(targetRxnIdx)/abs(mutSol_y(csRxnIdx));
bYield   = mutSol_y(bioRxnIdx)/(abs(mutSol_y(csRxnIdx))*csMW);
disp('Finding a minimal combination of targets displaying:')
disp([' - a production rate of: ' num2str(OptprodR) ' mmol/gDwh'])
disp([' - a production yield of: ' num2str(Optyield) ' mmol/mmol '])
disp([' - a biomass yield of: ' num2str(bYield) 'g biomass/g '])
disp(' ')
%sort targets by priority level and k_score
[levelCandidates,~] = sortrows(candidates,{'priority' 'k_scores'},{'ascend' 'ascend'});
counter   = 0;
remaining = table();
for i=1:height(levelCandidates)
    %reverse modifications
    gene   = levelCandidates.genes(i);
    action = levelCandidates.actions(i);
    enzIdx = levelCandidates.enz_pos(i);
    short  = levelCandidates.shortNames{i};
    mutF = 1;
    tempMutant = optMutant;
    %revert mutation
    if enzIdx>0
        saturationOpt         = mutSol_r(enzIdx)/(optMutant.lb(enzIdx)+-1e-15);
        tempMutant.lb(enzIdx) = 1.01*model.lb(enzIdx);
        tempMutant.ub(enzIdx) = 0.99*model.ub(enzIdx);
        if tempMutant.ub(enzIdx) <=tempMutant.lb(enzIdx)
           % tempMutant.lb(enzIdx) = 0.95*tempMutant.ub(enzIdx);
        end
        
    %for reactions without enzymatic reaction
    else
        enzUsage = 1e-12;
        if strcmpi(action,'KO')
            reversal = {'OE'};
        else
            if strcmpi(action,'KD')
                reversal = {'OE'};
            else
                reversal = {'KD'};
            end
        end
        modifications = {gene reversal mutF};
        tempMutant    = getMutantModel(optMutant,modifications,enzUsage);
        saturationOpt = NaN;
    end
    %Get max WT production rate
    mutsol_prod = solveECmodel(tempMutant,tempMutant,'pFBA','prot_pool_exchange');
    %Get max WT production yield
    mutsol_yield = solveECmodel(tempMutant,tempMutant,'pFBA',csRxn);
    mutprodR     = mutsol_prod(targetRxnIdx);
    mutyield     = mutsol_yield(targetRxnIdx)/abs(mutsol_yield(csRxnIdx));
    flag = true;
    if enzIdx>0 && numel(enzIdx)==1
        saturationM =  mutsol_yield(enzIdx)/(tempMutant.lb(enzIdx)+-1e-15);
        if (strcmpi(action,'OE') & saturationM <= (abs(model.lb(enzIdx))/levelCandidates.maxUsage(i)))
            flag = false;
        end
    else
        saturationM = NaN;
    end
    FC_y  = mutyield/Optyield;
    FC_p  = mutprodR/OptprodR;
    score = mean([FC_y,FC_p]);
    % Discard genes that don't affect the optimal phenotype
    % Discard OE targets that show low saturation after reversing the modification
    
    if isnan(score)
        score = 0;
    end
    
    if flag && ...
       (...
        (score<=1+tol) || ...
        (strcmpi(action,'OE') && saturationOpt >= 0.99) ...%%|| ((~strcmpi(action,'OE') && saturationOpt <= saturationM)) ...
       )
        
        remaining = [remaining;levelCandidates(i,:)];
        counter = counter+1;
        disp(['  Validated optimal target # ' num2str(counter) ': (' short '), ' action{1}])
    else
    	optMutant = tempMutant;
    end
end
remaining = sortrows(remaining,{'priority' 'k_scores'},{'ascend' 'descend'});
remaining = removevars(remaining,{'enz_pos' 'OE' 'minUsage' 'maxUsage' 'enzUsage' 'minUsageBio' 'maxUsageBio' 'enzUsageBio'});
end