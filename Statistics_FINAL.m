%% Statistics for Chapter 3: Hindlimb Proprioception

% Read in Data
D_all = readtable('Analyzed_Data_Final.csv');

% Only care about coordinated landings (exclude crash landings)
Data = D_all(D_all.Crash==0,:);

% Need to aggregate pseudoreplicate jump trials to individual averages for
% each condition
conditions = {'pre'; 'sham'; 'post'};
for cond = 1:size(conditions,1)
    Condition_Data{cond,1} = Data(strcmp(Data.condition,conditions{cond,1}),:);
    toads = unique(Condition_Data{cond,1}.toad);
    for indv = 1:size(toads,1)
        indv_Data{cond,indv} = Condition_Data{cond,1}(Condition_Data{cond,1}.toad==toads(indv,1),:);
        if cond == 1 && indv == 1
            D = [indv_Data{1,1}(1,[3,5]), varfun(@mean, indv_Data{1,1}(:,6:end))];
                else
            D = [D; [indv_Data{cond,indv}(1,[3,5]), varfun(@mean, indv_Data{cond,indv}(:, 6:end))]];
                end
    end
end

% Creating categorical variable for surgical condition
condNum = nan(size(D,1),1);
condNum(strcmp(D.condition,'pre')) = 1;
condNum(strcmp(D.condition,'sham')) = 2;
condNum(strcmp(D.condition,'post')) = 3;
D.condNum = categorical(condNum);
D.toad = categorical(D.toad); % define toad (individual) as categorical variable

%% Questions of Interest

% 1) Does elbow extension differ after nerve reinnervation? (extension rate, onset timing, & extension at
%    touchdown)
% 2) Is hindlimb/ankle extension affected by nerve reinnervation? (extension rate,
%    duration, & at touchdown)
% 3) Do activity of the Anconeus & Plataris change as a result of nerve
%    reinnervation? (timing of onset of plantaris prior to takeoff & anconeus prior to landing/after takeoff)
%    *Answered on Statistics_Q3.m file since using one-sample t-tests.

%% Statistical Rationale 
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

%% Statistical Analysis
modelStrings = { ' ~ 1 + (1|toad)' ... % model only containing toad as random effect
                 ' ~ 1 + condition + (1|toad)'}; % model with 1st order effects including condition as fixed effect

statsVars = D.Properties.VariableNames([3:13]); % [3:13] only compares variables of interest

logLkhdMatrix = nan(length(modelStrings)-1,length(statsVars));

for j = 1:length(statsVars) % loop through all dependent variables
    cVname = statsVars{j};
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

[bestModelAIC, bestModelIDx] = min(AICmatrix,[],1); % selecting model with lowest AIC

% Select the best model
modelToUse = mode(bestModelIDx); % select model that produced the most minmum AIC across all dependent variables 

% Retrieve stats from best model
for j = 1:length(statsVars)
    stats{1,j} = anova(compLME{modelToUse,j});
    rSquared(j) = compLME{modelToUse,j}.Rsquared;
    compLME{modelToUse,j};
    disp(statsVars{j});
    disp(stats{1,j})
    [b_rand{1,j},bnames{1,j},randStats{1,j}] = randomEffects(compLME{modelToUse,j});
end

%% Creating Table 1
for j = 1:size(statsVars,2)
    Variable = {statsVars{1,j}(6:end)};
    Var_idx = find(strcmp(D.Properties.VariableNames, statsVars{1,j}), 1);
    for c = 1:size(conditions,1)
        Condition = conditions(c,1);
        y = D{strcmp(D.condition,Condition{1,1}), Var_idx};
        Mean{1,j}(c,1) = mean(y, 'omitnan');
        SD{1,j}(c,1) = std(y, 'omitnan');
        CI_95prcnt{1,j}(c,1) = 1.96*SD{1,j}(c,1)/sqrt(size(rmmissing(y,1),1));
        F_stat = stats{1,j}{2,2};
        p_value = stats{1,j}{2,5};
        if j == 1 && c == 1
            Table_1 = table(Variable, Condition, Mean{1,j}(c,1), CI_95prcnt{1,j}(c,1), SD{1,j}(c,1), F_stat, p_value);
        elseif c == 1 % only care about variable name, F-statistic & p-value once for each variable
            Table_1 = [Table_1; table(Variable, Condition, Mean{1,j}(c,1), CI_95prcnt{1,j}(c,1), SD{1,j}(c,1), F_stat, p_value)];
        else % don't include redundant information
            Variable = {' '};
            F_stat = {' '};
            p_value = {' '};
            Table_1 = [Table_1; table(Variable, Condition, Mean{1,j}(c,1), CI_95prcnt{1,j}(c,1), SD{1,j}(c,1), F_stat, p_value)];
        end
    end
end


%% Visualizing Data
condColors = [0, 0.65, 1; % Light Blue
             0, 0.45, 0.75; % Medium Blue
             0, 0, 0.62000]; % Dark Blue

condCatVals = unique(condNum);
condition = condNum;
fontProperties = {'FontName','Arial','FontSize',14};

%% START FROM HERE!!!
for j = 1:length(statsVars)
    %Remove the effect of individual to run posthoc comparison of the fixed effect alone
    yFit_c = fitted(compLME{modelToUse,j},'Conditional',false); % Fit including fixed and random effects
    yFit_m = fitted(compLME{modelToUse,j},'Conditional',true); % Contribution from only fixed effects (marginal)
    y_ind = yFit_m - yFit_c; %Take the difference between the two above to get the effect of individual only. 
    
    %Calculate the corrected y-variable,  with the effect of individual removed
    cData_res = D.(statsVars{1,j}) - y_ind; 
    
    % Run an ANOVA model with 'anovan' to get the outputs necessary 
    % for the posthoc pairwise comparisons
   [pVals_ANOVAN{j},anovaTable{j},anovaStats{j}, anovaTerms{j}] = anovan(cData_res,...
       {D.condNum}...
        ,'random',[], 'continuous',[], 'model', 'interaction', 'varnames',...
        {'condNum'}, 'display', 'off');
    [pairwiseComps_c,MeansSE_c,~,groupNames_c] = multcompare(anovaStats{j},'ctype', 'bonferroni','Dimension',[1],'display','off');

    % Get pairwise mean differences and 95% confidence intervals for condition
    pairWisePvals{j} = [pairwiseComps_c(1:3,6)];
          
%% raincloud plot
    figH= figure;
    set(figH,'Color',[1 1 1]);
    hold on;
    cData = Data.(statsVars{1,j}(6:end));
    sortedData = cell(length(condCatVals),3);
    
    for k = 1:length(condCatVals)
        for m = 1:length(condCatVals)
            sortedData{k,m} = cData(strcmp(Data.condition,conditions{m,1}));
        end
    end
%   fh1{j} ORIGINALLY FIRST INPUT
    [ks_offsets(j,:)] = rm_raincloud_NO_IndvVariation(sortedData,condColors, 0, 'ks',CI_95prcnt{1,j},Mean{1,j});
%     [h,ks_offsets] = rm_raincloud(data, colours, plot_top_to_bottom, density_type, bandwidth)
%     for m = 1:size(fh1{i}.l,1)
%     fh1{i}.l(m,1).Color = [0.2 0.2 0.2];
%     fh1{i}.l(m,2).Color = [0.4 0.4 0.4]; 
%     end
%     for m = 1:size(fh1{i}.m,1)
%     fh1{i}.m(m,1).MarkerEdgeColor = [0.2 0.2 0.2];
%     fh1{i}.m(m,2).MarkerEdgeColor = [0.4 0.4 0.4];   
%     end
%     
    title(statsVars{1,j}(6:end));  
    %ylim([-0.5 6]);
    
        set(gca,'YTickLabel',fliplr({'pre' 'sham' 'post'}),'YTickLabelRotation',45)
        set(gca,'Box', 'off',fontProperties{:})
c_ax = gca;
set(c_ax,'XtickMode', 'manual','YTickMode','manual','Box','off','Color','none','FontName','Arial');

   set(figH,'units', 'centimeters');

%    width = 20;
%    height =  18.7500; %25
%    set(figH,'Position', [2 0.1 width height],'PaperPositionMode','auto')
   
   set(figH,'InvertHardcopy','off');
   figFilename = ['C:\Users\AJD44\Desktop\Chapter 3\Code\Figures\' statsVars{1,j}(6:end) '_RainCloud_' date '.eps'];
%    figure(figH)
%    saveas(figH,figFilename,'epsc')
% close(figH)
 end