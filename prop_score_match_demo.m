% restart
close all; clear; clc;

% load data
excel_filename = 'life_exp_data.xlsx';
all_data = readtable(excel_filename);

% compute probability of being in ctrl or treatment group based on covariates
% using logisitc regression 
x = [all_data.sex all_data.income];
treat = all_data.year;
[B,dev,stats] = mnrfit(x,treat);
prop = mnrval(B,x);
prop = prop(:,2);

% extract treatment and control propensity scores and outcome measures
treat_mask = (treat == 2);
ctrl_mask = ~treat_mask;
prop_treat = prop(treat_mask);
prop_ctrl = prop(ctrl_mask);
le_treat = all_data.life_exp(treat_mask);
le_ctrl = all_data.life_exp(ctrl_mask);

% compute raw (unadjusted) mean difference in life expectancy
raw_mean_le_treat = mean(le_treat);
raw_mean_le_ctrl = mean(le_ctrl);
raw_le_diff = raw_mean_le_treat-raw_mean_le_ctrl

% now do matching!
% using simple absolute deviation as the metric
matchedCtrlIdx = zeros(size(prop_treat)); % for each individual in treatment group, an index to the matched control
for treatIdx = 1:size(matchedCtrlIdx,1)
    diffs = abs(prop_ctrl-prop_treat(treatIdx));
    [~,minIdx] = min(diffs);
    matchedCtrlIdx(treatIdx) = minIdx;
end

% compute difference in mean life expectancy adjusted for covariates
adj_mean_le_treat = mean(le_treat);
adj_mean_le_ctrl = mean(le_ctrl(matchedCtrlIdx));
adj_le_diff = adj_mean_le_treat-adj_mean_le_ctrl