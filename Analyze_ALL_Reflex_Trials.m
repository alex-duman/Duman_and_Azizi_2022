% Analyzes all muscle stretch reflex trials & compile table of results for each trial 
Plot = 'N'; % set to 'Y' or 'Yes' to have all figures output, otherwise set to any other string or integer

% Read in a list of all trials
list = readtable('List_of_Reflex_Trials.xlsx');

% Calculates (and plots) variables for each trial & appends them to Table
rows = 6:size(list,1); % start from 6th file because first 5 are data from empty setup for calibration
for i = rows
    T = Analysis_and_Plot_Reflexes(char(list{i,1}),Plot);
    if i == min(rows)
        Table = T;
    else
        Table = [Table; T];
    end
end

writetable(Table, 'Analyzed_Data_Reflex.csv')
disp('Finished analyzing all reflex trials!')
disp('Results are saved within Analyzed_Data_Reflex.csv!')