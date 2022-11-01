% This script will determine and plot to check for a correlation between
% jump distance and elbow angle (extension) at touchdown.

% Read in the data
Data = readtable('Analyzed_Data_Final.csv');

% Breaking Data up into Pre, Sham and Post conditions
Conditions = unique(string(Data.condition));
Individuals = unique(Data.toad);
for i = 1:size(Conditions,1)
    D_cond{i,1} = Data(strcmp(string(Data.condition),Conditions(i,1)),:);

    for j = 1:size(Individuals,1)
        D{i,j} = D_cond{i,1}(D_cond{i,1}.toad==Individuals(j,1),:);
    end
end

Col = [1,0,0; 0,0,0; 0,0,1];
symbols = {'x'; 'd'; '*'; 'o'; '+'; 's'; 'v'; '<'; '>'; 'p'; '^'};
% Statistics and Plots
figure('Name','Relationship between Jump Distance & Elbow Extension')
for i = 1:size(Conditions,1)
    subplot(3,1,i)

    for j = 1:size(Individuals,1)
        data = rmmissing([D{i,j}.jump_dist, D{i,j}.ElbExt_td]);
        [R{i,j},P{i,j}] = corrcoef(data(:,1),data(:,2)); % Pearson's correlation
        r(i,j) = R{i,j}(1,2);
        p_corr(i,j) = P{i,j}(1,2);
        plot(data(:,1), data(:,2), symbols{j,1}, 'MarkerFaceColor',Col(i,:), 'MarkerEdgeColor', 'none');
        hold on;
        % determining regresssion line to plot relationship
        coeff{i,j} = polyfit(data(:,1), data(:,2), 1);
        X = [min(data(:,1)), max(data(:,1))];
        Y = coeff{i,j}(1,1)*X + coeff{i,j}(1,2);
        plot(X,Y,'-','Color',Col(i,:));
    end
    xlim([0, 0.25])
    ylim([20, 120])
    title(Conditions{i,1});
    ylabel('Elbow Extension (deg)');
    xlabel('Jump Distance (m)');
end

% Using a one-sample t-test to determine whether average correlation
% coefficient significantly differs from zero for each condition
for i = 1:size(Conditions,1)
    [~,p(i,1)] = ttest(rmmissing(r(i,:)));
end

% for i = 1:size(Conditions,1)
%     for j = 1:size(Individuals,1)
%         if coeff{i,j}(1,1) == 0 && coeff{i,j}(1,2) == 0 % no data to compute
%             slope(i,j) = nan;
%         else % calculated a meaningful regression
%             slope(i,j) = coeff{i,j}(1,1);
%         end
%     end
%     [~,p_slope(i,1)] = ttest(rmmissing(slope(i,:)));
% end