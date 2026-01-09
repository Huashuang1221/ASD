%% ============================================================
%  LMM for PageRank centrality (AUC) at left lateral OFC (ROI=26)
%  Model: y ~ 1 + Group + Age + IQ + (1|SiteID)
%  Output:
%    1) LMM_pagerank_OFC.mat                 (full LinearMixedModel object)
%    2) LMM_pagerank_OFC_results.mat         (results struct, tables)
%    3) LMM_pagerank_OFC_results.csv         (tidy table for sharing)
%    4) LMM_pagerank_OFC_workspace_snapshot.mat (snapshot for reproducibility)
%  MATLAB: R2018b+ (fitlme required)
%% ============================================================

clear; clc;

%% -------- 0) Set working directory (EDIT THIS) ---------------
workdir = 'D:\AAAResults_ASD_SMRI\wang_jinhui_results\jinhui_rsults_code';  
cd(workdir);

%% -------- 1) Load required data -----------------------------
% Node_para

if exist('Node_para.mat','file')
    load('Node_para.mat');  
elseif exist('Net_para.mat','file')
    load('Net_para.mat');   
else
    error('Cannot find Node_para.mat or Net_para.mat in %s', workdir);
end

% AGE & IQ covariates
if ~exist('AGE_IQ_424.mat','file')
    error('Cannot find AGE_IQ_424.mat in %s', workdir);
end
load('AGE_IQ_424.mat');  

% Site table
load('SITE.mat');      

% Group label (ASD vs TD)
Group = [];
if exist('ASD_ID.mat','file')
    tmp = load('ASD_ID.mat');            
    fn = fieldnames(tmp);
    ID = tmp.(fn{1});
    Group = repmat("TD",424,1);
    Group(ID) = "ASD";
elseif exist('ADOS_ID.mat','file')
    tmp = load('ADOS_ID.mat');
    fn = fieldnames(tmp);
    ID = tmp.(fn{1});
    Group = repmat("TD",424,1);
    Group(ID) = "ASD";
elseif exist('id.mat','file')
    tmp = load('id.mat');               
    fn = fieldnames(tmp);
    ID = tmp.(fn{1});
    Group = repmat("TD",424,1);
    Group(ID) = "ASD";
else
    warning('No ID file found. Using default grouping: 1:207 ASD, 208:424 TD.');
    Group = repmat("TD",424,1);
    Group(1:207) = "ASD";
end
Group = categorical(Group);

%% -------- 2) Prepare SiteID (categorical) -------------------

varnames = SITE.Properties.VariableNames;
siteVarName = varnames{1};  % default: first column
SiteID = categorical(SITE.(siteVarName));

% Quick checks
assert(numel(SiteID)==424, 'SiteID length is not 424.');
assert(size(AGE_IQ_424,1)==424 && size(AGE_IQ_424,2)>=2, 'AGE_IQ_424 must be 424x2 (Age, IQ).');

Age = AGE_IQ_424(:,1);
IQ  = AGE_IQ_424(:,2);

%% -------- 3) Extract PageRank (AUC) at ROI=26 ----------------
roi_OFC = 26;


if ~isfield(Node_para,'bin') || ~isfield(Node_para.bin,'aprank')
    error('Node_para.bin.aprank not found. Please check Node_para structure.');
end

y = Node_para.bin.aprank(:, roi_OFC);

% Basic QC
ok = isfinite(y) & isfinite(Age) & isfinite(IQ) & ~isundefined(SiteID) & ~isundefined(Group);
nSub = sum(ok);
fprintf('Valid subjects for LMM (pagerank): %d / 424\n', nSub);

%% -------- 4) Build table and fit LMM ------------------------
T = table(y(ok), Group(ok), Age(ok), IQ(ok), SiteID(ok), ...
    'VariableNames', {'y','Group','Age','IQ','SiteID'});

% Set TD as reference (optional but recommended)
T.Group = reordercats(T.Group, {'TD','ASD'});

formula = 'y ~ 1 + Group + Age + IQ + (1|SiteID)';
lme_pagerank = fitlme(T, formula, 'FitMethod','REML');

disp(lme_pagerank);
anovaTbl = anova(lme_pagerank, 'DFMethod','Residual');

%% -------- 5) Extract fixed effect of Group(ASD) --------------
coefTbl = lme_pagerank.Coefficients;


row = contains(coefTbl.Name, 'Group') & contains(coefTbl.Name, 'ASD');
if ~any(row)
    warning('Cannot find Group ASD coefficient by name. Printing coefficient table for manual check.');
    disp(coefTbl);
    row = strcmp(coefTbl.Name, 'Group_ASD'); % final attempt
end

GroupEffect = coefTbl(row, :);

%% -------- 6) Save results (MAT + CSV + snapshot) -------------
Results = struct();
Results.metric        = 'pagerank (AUC)';
Results.roi_index     = roi_OFC;
Results.formula       = formula;
Results.n_included    = nSub;
Results.site_var      = siteVarName;
Results.fixed_effects = coefTbl;
Results.group_effect  = GroupEffect;
Results.anova         = anova(lme_pagerank, 'DFMethod','Residual');
Results.random_effects = covarianceParameters(lme_pagerank);

% Tidy table for CSV export (only key terms)
out = table();
out.Metric   = string(Results.metric);
out.ROI      = Results.roi_index;
out.Term     = string(GroupEffect.Name);
out.Estimate = GroupEffect.Estimate;
out.SE       = GroupEffect.SE;
out.tStat    = GroupEffect.tStat;
out.DF       = GroupEffect.DF;
out.pValue   = GroupEffect.pValue;
out.CI_Lower = GroupEffect.Lower;
out.CI_Upper = GroupEffect.Upper;

% Filenames
f_model   = fullfile(workdir, 'LMM_pagerank_OFC.mat');
f_resmat  = fullfile(workdir, 'LMM_pagerank_OFC_results.mat');
f_csv     = fullfile(workdir, 'LMM_pagerank_OFC_results.csv');
f_snap    = fullfile(workdir, 'LMM_pagerank_OFC_workspace_snapshot.mat');

save(f_model, 'lme_pagerank');
save(f_resmat, 'Results', 'out', 'coefTbl', 'GroupEffect');
writetable(out, f_csv);

% snapshot (for reproducibility; optional but good practice)
save(f_snap);

fprintf('\nDONE (pagerank). Files saved:\n');
fprintf('  %s\n', f_model);
fprintf('  %s\n', f_resmat);
fprintf('  %s\n', f_csv);
fprintf('  %s\n', f_snap);
