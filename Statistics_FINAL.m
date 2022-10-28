%% Statistics for Duman & Azizi (2022)
% The authors would like to acknowledge and thank Dr. Monica Daley for
% providing initial code for performing LME Models.

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
% 1) Does elbow extension differ after nerve reinnervation? 
%    (extension rate, onset timing, & extension at touchdown)
% 2) Is hindlimb/ankle extension affected by nerve reinnervation? 
%    (extension rate, duration, & at touchdown)
% 3) Do activity of the Anconeus & Plataris change as a result of nerve
%    reinnervation? (timing of onset of plantaris prior to takeoff &
%    anconeus prior to landing/after takeoff)

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

%% Statistical Analysis (Questions 1 & 2)
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

%% Creating Table 1 (Table_1)
T1_VarNames = {'Variable', 'Condition', 'Mean', '95%_CI', 'SD', 'DF1', 'DF2', 'F_stat', 'p_value'};
for j = 1:size(statsVars,2)
    Variable = {statsVars{1,j}(6:end)};
    Var_idx = find(strcmp(D.Properties.VariableNames, statsVars{1,j}), 1);
    for c = 1:size(conditions,1)
        Condition = conditions(c,1);
        y = D{strcmp(D.condition,Condition{1,1}), Var_idx};
        Mean{1,j}(c,1) = mean(y, 'omitnan');
        SD{1,j}(c,1) = std(y, 'omitnan');
        CI_95prcnt{1,j}(c,1) = 1.96*SD{1,j}(c,1)/sqrt(size(rmmissing(y,1),1));
        DF1 = stats{1,j}{2,3}; % degrees of freedom for F-statistic numerator
        DF2 = stats{1,j}{2,4}; % degrees of freedom for F-statistic denominator
        F_stat = stats{1,j}{2,2};
        p_value = stats{1,j}{2,5};
        if j == 1 && c == 1
            Table_1 = table(Variable, Condition, Mean{1,j}(c,1), CI_95prcnt{1,j}(c,1), SD{1,j}(c,1), DF1, DF2, F_stat, p_value, 'VariableNames',T1_VarNames);
        elseif c == 1 % only care about variable name, F-statistic & p-value once for each variable
            Table_1 = [Table_1; table(Variable, Condition, Mean{1,j}(c,1), CI_95prcnt{1,j}(c,1), SD{1,j}(c,1), DF1, DF2, F_stat, p_value, 'VariableNames', T1_VarNames)];
        else % don't include redundant information
            Variable = {' '};
            DF1 = {' '};
            DF2 = {' '};
            F_stat = {' '};
            p_value = {' '};
            Table_1 = [Table_1; table(Variable, Condition, Mean{1,j}(c,1), CI_95prcnt{1,j}(c,1), SD{1,j}(c,1), DF1, DF2, F_stat, p_value, 'VariableNames', T1_VarNames)];
        end
    end
end

%% Visualizing Data (Questions 1 & 2)
condColors = [0, 0.65, 1; % Light Blue
             0, 0.45, 0.75; % Medium Blue
             0, 0, 0.62]; % Dark Blue

condCatVals = unique(condNum);
condition = condNum;
fontProperties = {'FontName','Times New Roman','FontSize',12};

for j = 1:length(statsVars)
    %Remove the effect of individual to run posthoc comparison of the fixed effect alone
    yFit_c = fitted(compLME{modelToUse,j},'Conditional',false); % Fit including fixed and random effects
    yFit_m = fitted(compLME{modelToUse,j},'Conditional',true); % Contribution from only fixed effects (marginal)
    y_ind = yFit_m - yFit_c; %Take the difference between the two above to get the effect of individual only. 
    
    %Calculate the corrected y-variable,  with the effect of individual removed
    cData_res = D.(statsVars{1,j}) - y_ind; 
    
    % Run an ANOVA model with 'anovan' to get the outputs necessary 
    % for the posthoc pairwise comparisons
   [pVals_ANOVAN{j},anovaTable{j},anovaStats{j}, anovaTerms{j}] = anovan(...
       cData_res,{D.condNum},'random',[], 'continuous',[], 'model',...
       'interaction', 'varnames',{'condNum'}, 'display', 'off');
    [pairwiseComps_c{1,j},MeansSE_c,~,groupNames_c{1,j}] = multcompare(anovaStats{j},...
        'ctype', 'bonferroni','Dimension',[1],'display','off');
          
% setting up raincloud plot
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

    Plot_RainClouds(sortedData,condColors, 12, 0, 'ks',CI_95prcnt{1,j},Mean{1,j});

    title(['(LME Model, F = ', num2str(round(stats{1,j}{2,2},2)),', p = ', num2str(round(stats{1,j}{2,5},3)), ')']); 
    xlabel(statsVars{1,j}(6:end)); % label for y-axis since axes are flipped 
    set(gca,'YTickLabel',fliplr({'pre' 'sham' 'post'}),'YTickLabelRotation',45)
    c_ax = gca;
    set(c_ax,'Box','on','Color','none','FontName','Times New Roman');
   
    set(figH,'InvertHardcopy','off');
    figFilename = ['C:\Users\AJD44\Desktop\Chapter 3\Code\Figures\' statsVars{1,j}(6:end) '_RainCloud_' date '.eps'];
%     saveas(figH,figFilename,'epsc')
%     close(figH)
end

%% Creating Table for Multiple Comparisons (Table_1_MultComp)
Condition_Comparisons = {'pre-sham'; 'pre-post'; 'sham-post'}; % determined from groupNames_c & pairwiseComps_c
for j = 1:size(statsVars,2)
    Variable = {statsVars{1,j}(6:end)};
    LME_p_value = stats{1,j}{2,5};
    MultComp_p_values = pairwiseComps_c{1,j}(1:3,6);

        if j == 1 % create new table
            Table_1_MultComp = table(Variable, LME_p_value, MultComp_p_values(1,1), MultComp_p_values(2,1), MultComp_p_values(3,1), ...
                'VariableNames', {'Variable', 'LME p-value', 'pre-sham p-val', 'pre-post p-val', 'sham-post p-val'});
        else % add to table
            Table_1_MultComp = [Table_1_MultComp; ...
                table(Variable, LME_p_value, MultComp_p_values(1,1), MultComp_p_values(2,1), MultComp_p_values(3,1), ...
                'VariableNames', {'Variable', 'LME p-value', 'pre-sham p-val', 'pre-post p-val', 'sham-post p-val'})];
        end
end

%% Statistical Analysis (Questions 3)
% Remove Pre-surgery and Sham trials
Data_EMG = Data(strcmp(Data.condition,'post')==1,[1:5,17:end]);
% Remove any rows of all NaNs for EMG data (data was not usable)
Data_EMG = Data_EMG(isnan(Data_EMG.Plant_act_dur_to)==0 | isnan(Data_EMG.Anc_ON_HLliftOff)==0 | isnan(Data_EMG.Anc_act_dur_air)==0,:);

% Need to aggregate pseudoreplicate jump trials to individual averages
DepVar = Data_EMG.Properties.VariableNames(6:8)';
for DV = 1:size(DepVar,1) % number of dependent variables for EMG
    DV_idx = find(strcmp(Data_EMG.Properties.VariableNames, DepVar{DV,1}), 1);
    DV_Data = Data_EMG(isnan(Data_EMG{:,DV_idx})==0, [5,DV_idx]); % only need individual and dependent variable of interest for comparison (only have EMG data for one condition, post-reinnervation)
    toads = unique(DV_Data.toad);
    for indv = 1:size(toads,1) % across all toads
        indv_Data{DV,indv} = DV_Data(DV_Data.toad==toads(indv,1),:);
        indv_mean = mean(indv_Data{DV,indv}{:,2},'omitnan');
        if indv == 1
            D_EMG{DV,1} = table(toads(indv,1), indv_mean, 'VariableNames', {'toad', ['mean_', DepVar{DV,1}]});
        else
            D_EMG{DV,1} = [D_EMG{DV,1};  table(toads(indv,1), indv_mean, 'VariableNames', {'toad', ['mean_', DepVar{DV,1}]})];
        end
    end
end

%% Average Values reported from literature in seconds
% Here, mu (µ) represents the average value(s) reported in the literature,
% but Matlab does not support the character µ for starting variable names.

% Duration of Plantaris activation during takeoff phase (Plant_act_dur_to)
mu_plant_dur_to = 0.130; % (Gillis & Biewener 2000; results on Plantaris)

% Initial Anconeus Activation during aerial phase relative to hindlimb liftoff (Anc_ON_touchdown)
mu_anc_on_to(1,1) = -0.048; % SD = 0.019, n_toads = 6, n_hops = 56 (Cox & Gillis 2020; Table 2)
mu_anc_on_to(1,2) = -0.036; % SD = 0.031, n_toads = 14, n_hops = 188 (Cox et al. 2018; Table 2)

% Duration of Anconeus activation during aerial phase (Anc_act_dur_air)
mu_anc_dur_air(1,1) = 0.088; % SD = 0.019, n_toads = 6, n_hops = 56 (Cox & Gillis 2020; Table 2)
mu_anc_dur_air(1,2) = 0.10; % SD = 0.032, n_toads = 14, n_hops = 188 (Cox et al. 2018; Table 2)

% Combine into single table mirroring order of Data table
Lit_vals = table(mu_plant_dur_to,mu_anc_on_to,mu_anc_dur_air);
Mean_EMG{1,1}(1,1) = mu_plant_dur_to;
Mean_EMG{1,2} = mu_anc_on_to';
Mean_EMG{1,3} = mu_anc_dur_air';

%% Calculate 95% Confidence Interval for Literature Values
% 95% CI = mean +/- 1.96*SD/sqrt(n), only need range or what is right of +/-

% Duration of Plantaris Activation during takeoff phase (Plant_act_dur_to)
plant_dur_CI = 0; % not enough information present in paper

% Initial Anconeus Activation during aerial phase relative to hindlimb liftoff (Anc_ON_touchdown)
anc_on_CI(1,1) = 1.96*0.019/sqrt(6); % SD = 0.019, n_toads = 6, n_hops = 56 (Cox & Gillis 2020)
anc_on_CI(1,2) = 1.96*0.031/sqrt(14); % SD = 0.031, n_toads = 14, n_hops = 188 (Cox et al. 2018)

% Duration of Anconeus activation during aerial phase (Anc_act_dur_air)
anc_dur_CI(1,1) = 1.96*0.019/sqrt(6); % SD = 0.019, n_toads = 6, n_hops = 56 (Cox & Gillis 2020)
anc_dur_CI(1,2) = 1.96*0.032/sqrt(14); % SD = 0.032 n_toads = 14, n_hops = 188 (Cox et al. 2018; Table 2)

% Combine into table mirroring Lit_vals
Lit_CIs = table(plant_dur_CI,anc_on_CI,anc_dur_CI);
CI{1,1}(1,1) = plant_dur_CI;
CI{1,2}(1,1) = anc_on_CI(1,1);
CI{1,2}(2,1) = anc_on_CI(1,2);
CI{1,3}(1,1) = anc_dur_CI(1,1);
CI{1,3}(2,1) = anc_dur_CI(1,2);


%% Combining Experimental Data with Literature Values
for j = 1:size(Lit_vals,2)
    if isnan(Lit_vals{1,j}) ~= 1 % if not a nan perform stats
        ci{1,j}(1:2,1) = nan;
        Mean_EMG{1,j}(size(Mean_EMG{1,j},1)+1, 1) = mean(D_EMG{j,1}{:,2}, 'omitnan');
        for k = 1:size(Lit_vals{1,j},2)
            [h{k,j},p{k,j},ci{k,j}] = ttest(D_EMG{j,1}{:,2},Lit_vals{1,j}(1,k));
            if k == 1 % only need to add CI for Experimental data once
                CI{1,j}(size(CI{1,j},1)+1, 1) = (ci{k,j}(2,1) - ci{k,j}(1,1))/2;
            end
        end
    else
        h{1,j} = nan;
        p{1,j} = nan;
        ci{1,j}(1:3,1) = nan;
    end
end

% Creating Table for EMG Statistical Results (Table_EMG)
TEMG_VarNames = {'Variable', 'Mean_post', '95%_CI_post', 'Mean_Ref1', '95%_CI_Ref1', 'p_val_Ref1', 'Mean_Ref2', '95%_CI_Ref2', 'p_val_Ref2'};
for j = 1:size(Lit_vals,2)
    Variable = DepVar(j,1);
    if j == 1 % create new table & we know we only have 1 reference (Ref1) & were unable to determine CI from paper
        Table_EMG = table(Variable, Mean_EMG{1,j}(end,1), CI{1,j}(end,1), Mean_EMG{1,j}(1,1), {' '}, p{1,j}, {' '}, {' '}, {' '}, 'VariableNames',TEMG_VarNames);
    else % add to table
        Table_EMG = [Table_EMG; table(Variable, Mean_EMG{1,j}(end,1), CI{1,j}(end,1), Mean_EMG{1,j}(1,1), CI{1,j}(1,1), p{1,j}, Mean_EMG{1,j}(2,1), CI{1,j}(2,1), p{2,j},'VariableNames',TEMG_VarNames)];
    end
end

%% Visualizing Data
fontProperties = {'FontName','Times New Roman','FontSize',12};
condColors_EMG = [1, 0.95, 0; % Yellow
                  0.75, 0.46, 0.16; % Burnt Orange
                  0, 0, 0.62]; % Dark Blue (for post-reinnervation)

condNum = zeros(size(Data_EMG,1),1) + 1; % always 'post' condition so only need 1 as condition number
Data_EMG.condNum = categorical(condNum);
condCatVals = unique(condNum);

for j = 1:size(Lit_vals,2)
    cVname = Data_EMG.Properties.VariableNames{j+5};
    figH= figure;
    set(figH,'Color',[1 1 1]);
    hold on;
    cData = Data_EMG.(cVname);
    sortedData = cell(length(condCatVals),(size(Lit_vals{1,j},2)+1));
    
    for k = 1:(size(Lit_vals{1,j},2)+1)
        sortedData{k,1} = Lit_vals{1,j}(1,1);
        if size(Lit_vals{1,j},2)==2
            sortedData{k,2} = Lit_vals{1,j}(1,2);
            sortedData{k,3} = rmmissing(cData);
        else
            sortedData{k,2} = rmmissing(cData);
        end
    end

    if size(Lit_vals{1,j},2)==2 % have two reference values from literature
        Plot_RainClouds(sortedData,condColors_EMG, 12, 0, 'ks', CI{1,j},Mean_EMG{1,j})
        set(gca,'YTickLabel',fliplr({'Ref1' 'Ref2' 'post'}),'YTickLabelRotation',45)
        title({['(Ref1, t-test, p =', num2str(round(p{1,j},3)), ')']; ['(Ref2, t-test, p =', num2str(round(p{2,j},3)), ')']})
    else % only have one reference value from literature
        Plot_RainClouds(sortedData,condColors_EMG(2:3,:), 12, 0, 'ks', CI{1,j},Mean_EMG{1,j})
        set(gca,'YTickLabel',fliplr({'Ref1' 'post'}),'YTickLabelRotation',45)
        title(['(t-test, p = ', num2str(round(p{1,j},3)), ')'])
    end
    xlabel(cVname); % label for y-axis since axes are flipped
    c_ax = gca;
    set(c_ax,'Box','on','Color','none','FontName','Times New Roman');
   
    set(figH,'InvertHardcopy','off');
    figFilename = ['C:\Users\AJD44\Desktop\Chapter 3\Code\Figures\' cVname '_RainCloud_' date '.eps'];
    figure(figH)
%     saveas(figH,figFilename,'epsc')
 end

%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Plot_RainClouds(data, colours, sz, plot_top_to_bottom, density_type, errorBars,Means, bandwidth)
% This code is originally modified from RainCloudPlots/RainCloudPlots
% GitHub repository found at: 
% https://github.com/RainCloudPlots/RainCloudPlots

% Use like: Plot_RainClouds(data, colours, plot_top_to_bottom, density_type, errorBars, Means, bandwidth)

% INPUTS
% data - an M x N cell array of N data series and M measurements
% colours - an N x 3 array defining the colour to plot each data series
% sz - scalar of marker size (e.g. 12)
% plot_top_to_bottom - Default plots left-to-right, set to 1 to rotate.
% density_type - 'ks' (default) or 'RASH'. 'ks' uses matlab's inbuilt 
%                'ksdensity' to determine the shape of the rainclouds. 
%                'RASH' will use the 'rst_RASH' method from Cyril Pernet's 
%                Robust Stats toolbox, if that function is on your matlab path.
% errorBars - vector for error bar values (e.g. standard deviation or confidence
%             intervals, etc.)
% Means - are the mean values calculated for each group as a vector
% bandwidth - (optional) if density_type == 'ks', determines bandwidth of density estimate


% check dimensions of data
[nper, nseries] = size(data);

% make sure we have enough colours
assert(all(size(colours) == [nseries 3]), 'number of colors does not match number of plot series');

if nargin < 6
    plot_top_to_bottom = 0; % left-to-right plotting by default
end

if nargin < 7
    density_type = 'ks';
end

if nargin < 8
    bandwidth = [];
end

%% Calculate properties of density plots

% Probably okay to hard-code this as it just determines the granularity of
% the density estimate

nbins = repmat(100, nper, nseries);

% calculate kernel densities
for i = 1:nper
    for j = 1:nseries
        % TO-DO: Switch here to use alternative methods to estimate density (e.g. RASH)
        switch density_type
            case 'ks'
                [ks{i,j}, x{i,j}] = ksdensity(data{i,j}, 'NumPoints', nbins(i,j), 'bandwidth', bandwidth);
            case 'rash'
                % check for rst_RASH function (from Robust stats toolbox) in path, fail if not found 
                assert(exist('rst_RASH','file') == 2, 'Could not compute density using RASH method. Do you have the Robust Stats toolbox on your path?');
                
                [x{i,j}, ks{i,j}] = rst_RASH(data{i,j});
                % override default 'nbins' as rst_RASH determines number of bins
                nbins(i,j) = size(ks{i,j},2);
        end
        
        % Define the faces to connect each adjacent f(x) and the corresponding points at y = 0.
        q{i,j}       = (1:nbins(i,j)-1)';
        faces{i,j}   = [q{i,j}, q{i,j} + 1, q{i,j} + nbins(i,j) + 1, q{i,j} + nbins(i,j)];
    end
end

% keyboard

% determine spacing between plots
spacing     = 2 * mean(mean(cellfun(@max,ks)));
ks_offsets  = [0:nper-1] .* spacing;
% flip so first plot in series is plotted on the *top*
ks_offsets  = fliplr(ks_offsets);

% calculate patch vertices from kernel density
for i = 1:nper
    for j = 1:nseries
        verts{i,j} = [x{i,j}', ks{i,j}' + ks_offsets(i); x{i,j}', ones(nbins(i,j),1) * ks_offsets(i)];
        verts{i,j} = [x{i,j}', ks{i,j}' + ks_offsets(i); x{i,j}', ones(nbins(i,j),1) * ks_offsets(i)];
    end
end


%% jitter for the raindrops

jitwidth = spacing / 8;

for i = 1:nper
    for j = 1:nseries
        jit{i,j} = jitwidth + rand(1, length(data{i,j})) * jitwidth;
    end
end

%% plot
% note - we *could* plot everything here in one big loop, but then
% different figure parts would overlay each other in a silly way.

hold on

% patches & raindrops
for i = 1:nper
    for j = 1:nseries
        if j == i % only plot data for specific condtion
            if size(data{i,j},1)>1
                % plot patches
%               plot(verts{i,j}(:,1),verts{i,j}(:,2),'-','Color',colours(j,:)); % use if only want outline of cloud distribution
                patch([verts{i,j}(1:100,1), fliplr(verts{i,j}(1:100,1))], [verts{i,j}(101:200,2), fliplr(verts{i,j}(1:100,2))], colours(j,:),'FaceAlpha',0.5, 'EdgeColor','none')

                % scatter plot of raindrops
                scatter(data{i,j}, -jit{i,j} + ks_offsets(i), sz, 'MarkerFaceColor', colours(j,:), 'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.5);
            end
        end
    end
end

% plot mean dots
for i = 1:nper
    for j = 1:nseries
        if i==j
            plot(Means(i,1),ks_offsets(i),'o','MarkerFaceColor', colours(j,:), 'MarkerEdgeColor', [0 0 0], 'MarkerSize', sz*0.6);
            % add error bars to mean point (e.g. standard deviation, 95% CI, etc.)
            er = errorbar(Means(i,1),ks_offsets(i), errorBars(i,1),'horizontal','CapSize',0,'LineWidth',2); er.Color = [0,0,0]; er.LineStyle = 'none';
            % Re-plot mean so they appear on top of error bars
            plot(Means(i,1),ks_offsets(i),'o','MarkerFaceColor', colours(j,:), 'MarkerEdgeColor', [0 0 0], 'MarkerSize', sz*0.6);
        end
    end
end

%% clear up axis labels

% 'YTick', likes values that *increase* as you go up the Y-axis, but we plot the first
% raincloud at the top. So flip the vector around
set(gca, 'YTick', fliplr(ks_offsets));
set(gca, 'YTickLabel', nper:-1:1);
ylim([(min(ks_offsets)-(spacing/2)) (max(ks_offsets)+spacing)]);

%% determine plot rotation
% default option is left-to-right
% plot_top_to_bottom can be set to 1 
% NOTE: Because it's easier, we actually prepare everything plotted
% top-to-bottom, then - by default - we rotate it here. That's why the
% logical is constructed the way it is.

% rotate and flip
if ~plot_top_to_bottom
    view([90 -90]);
    axis ij
end
end