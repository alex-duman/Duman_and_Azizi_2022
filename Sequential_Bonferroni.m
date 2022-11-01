LME_stats = readtable('Compiled_Statistics.xlsx','Sheet',1);
Ttest_stats = readtable('Compiled_Statistics.xlsx','Sheet',2);

Variables_LME = unique(LME_stats.Variable);
row = 1;
for i = 2:size(Variables_LME,1) % don't need first variable since this is all the blanks
    if row == 1 % start new table
        Seq_Bon = table({LME_stats.Variable{strcmp(LME_stats.Variable,Variables_LME{i,1})}(2:end-1)}, str2num(LME_stats.p_value{strcmp(LME_stats.Variable,Variables_LME{i,1})}),'VariableNames', {'Variable', 'p_value'});
        row = row+1;
    else % add to table
        Seq_Bon = [Seq_Bon; table({LME_stats.Variable{strcmp(LME_stats.Variable,Variables_LME{i,1})}(2:end-1)}, str2num(LME_stats.p_value{strcmp(LME_stats.Variable,Variables_LME{i,1})}), ...
            'VariableNames', {'Variable', 'p_value'})];
    end
end

Variables_Ttest = unique(Ttest_stats.Variable);
for i = 1:size(Variables_Ttest,1)
    Seq_Bon = [Seq_Bon; table({Ttest_stats.Variable{strcmp(Ttest_stats.Variable,Variables_Ttest{i,1})}(2:end-1)}, Ttest_stats.p_val_Ref1(strcmp(Ttest_stats.Variable,Variables_Ttest{i,1})), 'VariableNames', {'Variable', 'p_value'})];
    if ~isnan(Ttest_stats.p_val_Ref3(strcmp(Ttest_stats.Variable,Variables_Ttest(i,1)))) % if not missing then add to table
        Seq_Bon = [Seq_Bon; table({Ttest_stats.Variable{strcmp(Ttest_stats.Variable,Variables_Ttest{i,1})}(2:end-1)}, Ttest_stats.p_val_Ref3(strcmp(Ttest_stats.Variable,Variables_Ttest{i,1})), 'VariableNames', {'Variable', 'p_value'})];
    end
end

Seq_Bon_ranked = sortrows(Seq_Bon,"p_value");

for i = 1:size(Seq_Bon_ranked,1)
    Threshold(i,1) = 0.05/(size(Seq_Bon_ranked,1) - i + 1);
    if Seq_Bon_ranked.p_value(i) < Threshold(i,1)
        Significant{i,1} = 'True';
    else
        Significant{i,1} = 'False';
    end
end
% Set any threshold value
Seq_Bon_ranked = [Seq_Bon_ranked, table(Threshold, Significant)];