%% =========================================================
%  LMM: Betweenness centrality (left lateral OFC, ROI=26)
%  Random intercept: Site
%  Fixed effects: Group (ASD vs TD) + Age + IQ
%  MATLAB: R2018b
%% =========================================================

clear; clc;

%% 0) Path (official GRETNA root)
restoredefaultpath;
rehash toolboxcache;
addpath(genpath('F:\GRETNA_full'));
savepath;
rehash;

%% 1) Load
cd('D:\huashuang');   
load('Node_para.mat');      
load('AGE_IQ_424.mat');   
load('SITE.mat'); 

%% 2) Build predictors
nSub = size(AGE_IQ_424,1);
assert(nSub==424, 'Expect 424 subjects, got %d.', nSub);

Group = categorical([repmat({'ASD'},207,1); repmat({'TD'},217,1)]);
Age   = AGE_IQ_424(:,1);
IQ    = AGE_IQ_424(:,2);

siteVarName = SITE.Properties.VariableNames{1};
SiteID = categorical(SITE.(siteVarName));

disp(summary(SiteID)); 

%% 3) Outcome: betweenness AUC (abw), ROI=26
roi_OFC = 26;
assert(isfield(Node_para,'bin') && isfield(Node_para.bin,'abw'), 'Node_para.bin.abw not found.');
y = Node_para.bin.abw(:, roi_OFC);

% remove missing
ok = ~isnan(y) & ~isnan(Age) & ~isnan(IQ) & ~isundefined(SiteID) & ~isundefined(Group);
T = table(y(ok), Group(ok), Age(ok), IQ(ok), SiteID(ok), ...
          'VariableNames', {'y','Group','Age','IQ','SiteID'});

T.Group = reordercats(T.Group, {'TD','ASD'});  % TD as reference

%% 4) Fit LMM
formula = 'y ~ 1 + Group + Age + IQ + (1|SiteID)';
lme = fitlme(T, formula, 'FitMethod','REML');

disp(lme);
anovaTbl = anova(lme, 'DFMethod','Residual');
disp(anovaTbl);

coef = lme.Coefficients;
row = strcmp(coef.Name, 'Group_ASD');
if ~any(row)
    warning('Cannot find Group_ASD row. Names are:');
    disp(coef.Name);
end

%% 5) Save outputs
lme_betweenness = lme;
save('LMM_betweenness_OFC.mat', 'lme_betweenness');

Bet_LMM_Result = struct();
Bet_LMM_Result.Metric   = 'Betweenness centrality (AUC, Node_para.bin.abw)';
Bet_LMM_Result.ROI      = roi_OFC;
Bet_LMM_Result.Formula  = formula;
Bet_LMM_Result.N        = height(T);

if any(row)
    Bet_LMM_Result.Effect    = 'ASD vs TD (Group_ASD)';
    Bet_LMM_Result.Estimate  = coef.Estimate(row);
    Bet_LMM_Result.SE        = coef.SE(row);
    Bet_LMM_Result.tStat     = coef.tStat(row);
    Bet_LMM_Result.DF        = coef.DF(row);
    Bet_LMM_Result.pValue    = coef.pValue(row);
    Bet_LMM_Result.CI_Lower  = coef.Lower(row);
    Bet_LMM_Result.CI_Upper  = coef.Upper(row);
end

save('LMM_betweenness_OFC_results.mat', 'Bet_LMM_Result', 'coef', 'anovaTbl');
BwTable = struct2table(Bet_LMM_Result);
writetable(BwTable, 'LMM_betweenness_OFC_results.csv');

save('LMM_betweenness_OFC_workspace_snapshot.mat', 'T', 'roi_OFC');

fprintf('\nDONE (betweenness). Files saved:\n');
fprintf('  LMM_betweenness_OFC.mat\n');
fprintf('  LMM_betweenness_OFC_results.mat\n');
fprintf('  LMM_betweenness_OFC_results.csv\n');
fprintf('  LMM_betweenness_OFC_workspace_snapshot.mat\n');
