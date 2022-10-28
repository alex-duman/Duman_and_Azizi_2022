%% Statistics for Duman and Azizi (2022) Reflex Comparison

% Read in Data
Data = readtable('Analyzed_Data_Reflex.csv');

condNum = nan(size(Data,1),1);
condNum(strcmp(Data.Condition,'Pre-op')) = 1;
condNum(strcmp(Data.Condition,'Sham')) = 2;
condNum(strcmp(Data.Condition,'1wk Post-op')) = 3;
condNum(strcmp(Data.Condition,'3mo Post-op')) = 4;
condNum(strcmp(Data.Condition,'6mo Post-op')) = 5;
condNum(strcmp(Data.Condition,'Post-mortem')) = 6;
Data.condNum = categorical(condNum);

% Remove any trials that resulted in Ratios of 0 or Inf (either recorded no
% activity during flexion or extension phase which will affect individual averages)
for row = 1:size(Data,1)
    if Data{row,'Ratio_I_ext_flex'} == 0 || Data{row,'Ratio_I_ext_flex'} == Inf
        Data{row, 'Ratio_I_ext_flex'} = NaN;
    end
    if Data{row,'Ratio_t_on_ext_flex'} == 0 || Data{row,'Ratio_t_on_ext_flex'} == Inf
        Data{row, 'Ratio_t_on_ext_flex'} = NaN;
    end
end

% Need to aggregate data by getting the average ratio for each individual
% at every timepoint to remove pseudoreplication from dataset
conditions = {'Pre-op'; 'Sham'; '1wk Post-op'; '3mo Post-op'; '6mo Post-op'; 'Post-mortem'};
for cond = 1:size(conditions,1)
    Condition_Data{cond,1} = Data(strcmp(Data.Condition,conditions{cond,1}),:);
    toads = unique(Condition_Data{cond,1}.Individual);
    for indv = 1:size(toads,1)
        indv_Data{cond,indv} = Condition_Data{cond,1}(Condition_Data{cond,1}.Individual==toads(indv,1),:);
        if cond == 1 && indv == 1
            D = [indv_Data{1,1}(1,1:2), table(mean(indv_Data{1,1}.Ratio_I_ext_flex, 'omitnan'), mean(indv_Data{1,1}.Ratio_t_on_ext_flex, 'omitnan'), 'VariableNames',{'Ratio_I_ext_flex', 'Ratio_t_on_ext_flex'})];
        else
            D = [D; [indv_Data{cond,indv}(1,1:2), table(mean(indv_Data{cond,indv}.Ratio_I_ext_flex, 'omitnan'), mean(indv_Data{cond,indv}.Ratio_t_on_ext_flex, 'omitnan'), 'VariableNames',{'Ratio_I_ext_flex', 'Ratio_t_on_ext_flex'})]];
        end
    end
end

%% Questions of Interest

% 1) Does plantaris activity significantly differ before and after nerve
% reinnervation?
%       A) Does the ratio of EMG Intensity during Extension/Flexion change?
%       B) Does the ratio of Active Time during Extension/Flexion change?

%% Statistical Rationale (Questions 1 & 2)
% Linear mixed effects (LME) models 'fitlme' are better than ANOVAs in this
% instance since they deal better with unbalanced datasets (some toads 
% didn't participate in sham or post-surgery condition).

% Our models treat 'toad' as a random effect and the 'Condition' (pre, sham 
% and post-surgery) as fixed effect.
 
% We test the simplest model first, with no fixed effect and only the 
% random effect of individual variation being present (y ~ (1|toad)). We
% can then reference this model to test and compare to when we include the
% fixed effect of treatment condition (y ~ (1|toad) + condition). We
% compare the AIC and/or R-squared of the two models. We then look for
% the model that produces the lowest AIC and/or highest R-squared across
% the most dependent variables we investigate.

%% Statistical Analysis (Checking Muscle Stretch Reflex)
modelStrings = { ' ~ 1 + (1|Individual)' ... % model only containing toad as random effect
                 ' ~ 1 + Condition + (1|Individual)'}; % model with only 1st order effects

Cols = 3:4; % columns of the variables interested in
ReflexVars = D.Properties.VariableNames(Cols); % only compares variables of interest

logLkhdMatrix = nan(length(modelStrings)-1,length(ReflexVars));

for j = 1:length(ReflexVars) % loop through all dependent variables
    cVname = ReflexVars{j};
    for k = 1:length(modelStrings) % loop through all models
        cModelStr = modelStrings{k};
        compLME{k,j} = fitlme(D,[cVname cModelStr]);
        AICmatrix(k,j) = compLME{k,j}.ModelCriterion{1, 1};
        if k > 1
            t_loglklyhood = compare(compLME{k-1,j}, compLME{k,j});
            logLkhdMatrix(k-1,j) = t_loglklyhood.pValue(2);   
        end
    end
end

[bestModelAIC,bestModelIDx] = min(AICmatrix,[],1); %find model with lowest AIC

%% Get stats for the best model
modelToUse = mode(bestModelIDx); %find model that modal lowest AIC across variables 

for j = 1:length(ReflexVars)
    stats{1,j} = anova(compLME{modelToUse,j});
    rSquared(j) = compLME{modelToUse,j}.Rsquared;
    compLME{modelToUse,j} %#ok<NOPTS>
    disp(ReflexVars{j});
    disp(stats{1,j})
    [b_rand{1,j},bnames{1,j},randStats{1,j}] = randomEffects(compLME{modelToUse,j});
end

%% Creating Table of Statistical Results (Table_Reflex)
T_VarNames = {'Variable', 'Condition', 'Mean', '95%_CI', 'SD', 'DF1', 'DF2', 'F_stat', 'p_value'};
for j = 1:size(ReflexVars,2)
    Variable = {ReflexVars{1,j}(6:end)};
    Var_idx = find(strcmp(D.Properties.VariableNames, ReflexVars{1,j}), 1);
    for c = 1:size(conditions,1)
        Condition = conditions(c,1);
        y = D{strcmp(D.Condition,Condition{1,1}), Var_idx};
        Mean{1,j}(c,1) = mean(y, 'omitnan');
        SD{1,j}(c,1) = std(y, 'omitnan');
        CI_95prcnt{1,j}(c,1) = 1.96*SD{1,j}(c,1)/sqrt(size(rmmissing(y,1),1));
        DF1 = stats{1,j}{2,3}; % degrees of freedom for F-statistic numerator
        DF2 = stats{1,j}{2,4}; % degrees of freedom for F-statistic denominator
        F_stat = stats{1,j}{2,2};
        p_value = stats{1,j}{2,5};
        if j == 1 && c == 1
            Table_Reflex = table(Variable, Condition, Mean{1,j}(c,1), CI_95prcnt{1,j}(c,1), SD{1,j}(c,1), DF1, DF2, F_stat, p_value, 'VariableNames',T_VarNames);
        elseif c == 1 % only care about variable name, F-statistic & p-value once for each variable
            Table_Reflex = [Table_Reflex; table(Variable, Condition, Mean{1,j}(c,1), CI_95prcnt{1,j}(c,1), SD{1,j}(c,1), DF1, DF2, F_stat, p_value, 'VariableNames', T_VarNames)];
        else % don't include redundant information
            Variable = {' '};
            DF1 = {' '};
            DF2 = {' '};
            F_stat = {' '};
            p_value = {' '};
            Table_Reflex = [Table_Reflex; table(Variable, Condition, Mean{1,j}(c,1), CI_95prcnt{1,j}(c,1), SD{1,j}(c,1), DF1, DF2, F_stat, p_value, 'VariableNames', T_VarNames)];
        end
    end
end


%% Visualize Reflex Data
condColors = [0, 0.65, 1; % Light Blue (pre-op)
             0, 0.45, 0.75; % Medium Blue (sham)
             0, 0.40, 0.7; % blue (1wk post-op)
             0.04, 0.30, 0.63; % blue (3mo post-op)
             0, 0, 0.62000; % Dark Blue (6mo post-op)
             0.8, 0.8, 0.8]; % gray (post-mortem)

condCatVals = unique(condNum);
condition = condNum;

fontProperties = {'FontName','Times New Roman','FontSize',12};
x_labels = {'Pre-op', 'Sham (1wk)', '1wk Post-op', '3mo Post-op', '6mo Post-op', 'Post-mortem'};
for j = 1:length(ReflexVars)
    cVname = ReflexVars{j};   
    %Get the random effect results
    %Remove the effect of individual to run posthoc comparison of the fixed effect alone
    yFit_c = fitted(compLME{modelToUse,j},'Conditional',false); % Fit including fixed and random effects
    yFit_m = fitted(compLME{modelToUse,j},'Conditional',true); % Contribution from only fixed effects (marginal)
    y_ind = yFit_m - yFit_c; % Take the difference between the two above to get the effect of individual only. 
    
    %Calculate the corrected y-variable,  with the effect of individual removed
    cData_res = D.(cVname) - y_ind; 
    
    % Run an ANOVA model with 'anovan' to get the outputs necessary 
    % for the posthoc pairwise comparisons
   [pVals_ANOVAN{j},anovaTable{j},anovaStats{j}, anovaTerms{j}] = anovan(cData_res,...
       {D.Condition}...
        ,'random',[], 'continuous',[], 'model', 'interaction', 'varnames',...
        {'condNum'}, 'display', 'off');
    [pairwiseComps_c{j},MeansSE_c,~,groupNames_c] = multcompare(anovaStats{j},'ctype', 'bonferroni','Dimension',[1],'display','off');

    % Plot the means in a bar plot with 95% Confidence Intervals as error bars
    mean_values = Table_Reflex.Mean([1:size(conditions,1)]+size(conditions,1)*(j-1));
    CI_values = Table_Reflex.("95%_CI")([1:size(conditions,1)]+size(conditions,1)*(j-1));
    fh(j) = figure;
    set(gcf,'Color', [1 1 1]);
    b = bar(1:size(mean_values,1),mean_values, 'FaceColor', 'Flat');
    hold on;      
    er = errorbar(1:size(mean_values,1), mean_values, CI_values); er.Color = [0,0,0]; er.LineStyle = 'none';
    for i = 1:size(mean_values,1)
        b.CData(i,:) = condColors(i,:); % set color of each bar
    end
    ylabel(cVname);
    xlabel('Condition')
    title(['(LME Model, F = ', num2str(round(stats{1,j}{2,2},2)),', p = ', num2str(round(stats{1,j}{2,5},3)), ')'])
    set(gca,'XTickLabel',x_labels,'Box', 'on','XTickLabelRotation',45,fontProperties{:});
    
    set(fh(j),'InvertHardcopy','off');   
    figure(fh(j))
    hold off;
    filename = ['C:\Users\AJD44\Desktop\Chapter 3\Code\Figures\' cVname '_MeanDiff_' date '.eps'];
%     saveas(fh(j),filename,'epsc')
%     close(fh(j));
end

%% Creating Table for Multiple Comparisons (Table_1_MultComp)
Condition_Comparisons = {'pre-sham'; 'pre-post'; 'sham-post'}; % determined from groupNames_c & pairwiseComps_c
for j = 1:size(ReflexVars,2)
    for i = 1:size(pairwiseComps_c{1,j},1)
        Condition_Comparisons{i,1} = [groupNames_c{pairwiseComps_c{1,j}(i,1)}(9:end), ' - ', groupNames_c{pairwiseComps_c{1,j}(i,2)}(9:end)];
    end
    Variable = {ReflexVars{1,j}};
    LME_p_value = stats{1,j}{2,5};
    MultComp_p_values = pairwiseComps_c{1,j}(1:size(pairwiseComps_c{1,j},1),6);
    table_array = [LME_p_value, MultComp_p_values'];
    Temp_table = table(Variable, table_array);
        if j == 1 % create new table
            Table_Reflex_MultComp = splitvars(Temp_table, 'table_array', 'NewVariableNames', {'LME p-value', Condition_Comparisons{:,1}});
        else % add to table
            Table_Reflex_MultComp = [Table_Reflex_MultComp; ...
                splitvars(Temp_table, 'table_array', ...
                'NewVariableNames', {'LME p-value', Condition_Comparisons{:,1}})];
        end
end