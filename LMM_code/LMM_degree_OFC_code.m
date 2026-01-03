%% =========================================================
%  LMM: Degree centrality (left lateral OFC, ROI=26)
%  Random intercept: Site
%  Fixed effects: Group (ASD vs TD) + Age + IQ

%% =========================================================

clear; clc;


restoredefaultpath;
rehash toolboxcache;

addpath(genpath('F:\GRETNA_full'));
savepath;
rehash;

%% -------------------------
% 1) 载入数据
%    - Node_para: ComBat后网络参数
%    - AGE_IQ_424: 协变量（424×2）
%    - SITE: 站点信息（424×1 table，从excel导入的）
% --------------------------

cd('D:\huashuang');  

load('Net_para.mat');         
load('Node_para.mat');        % 必须：包含 Node_para
load('AGE_IQ_424.mat');       % 必须：AGE_IQ_424 (424×2)
SITE = readtable('SITE.xlsx');

%% -------------------------
% 2) 基本检查
% --------------------------
assert(exist('Node_para','var')==1, 'Node_para not found.');
assert(exist('AGE_IQ_424','var')==1, 'AGE_IQ_424 not found.');
assert(exist('SITE','var')==1, 'SITE table not found.');

nSub = size(AGE_IQ_424,1);
assert(nSub==424, 'Expect 424 subjects, but got %d.', nSub);

% 检查 Node_para.bin.adeg 尺寸
assert(isfield(Node_para,'bin') && isfield(Node_para.bin,'adeg'), 'Node_para.bin.adeg not found.');
assert(size(Node_para.bin.adeg,1)==nSub, 'Node_para.bin.adeg rows mismatch.');

%% -------------------------
% 3) 构造变量：Group / Age / IQ / SiteID
% --------------------------
% Group: 1-207 ASD, 208-424 TD
Group = categorical([repmat({'ASD'},207,1); repmat({'TD'},217,1)]);

% Age & IQ
Age = AGE_IQ_424(:,1);
IQ  = AGE_IQ_424(:,2);

% SiteID: 从 SITE table 取“第一列”作为站点标签（避免列名不一致问题）
siteVarName = SITE.Properties.VariableNames{1};  % 例如 'GU'
SiteID = categorical(SITE.(siteVarName));

% 快速检查
fprintf('N subjects = %d\n', nSub);
disp(summary(SiteID));   % 每个site人数

%% -------------------------
% 4) 提取目标指标：Degree centrality (OFC ROI = 26)
% --------------------------
roi_OFC = 26;  % 你之前 OFC ROI=26
y = Node_para.bin.adeg(:, roi_OFC);

% 检查缺失
if any(isnan(y))
    warning('y contains NaN. Removing rows with missing values.');
    ok = ~isnan(y) & ~isnan(Age) & ~isnan(IQ) & ~isundefined(SiteID) & ~isundefined(Group);
    y = y(ok);
    Age = Age(ok);
    IQ = IQ(ok);
    SiteID = SiteID(ok);
    Group = Group(ok);
end

%% -------------------------
% 5) 组装数据表并拟合 LMM
% --------------------------
T = table(y, Group, Age, IQ, SiteID);

% 建议固定参考水平为TD（更直观：Group_ASD 表示 ASD 相对 TD 的差异）
T.Group = reordercats(T.Group, {'TD','ASD'});

% LMM公式：随机截距 SiteID
formula = 'y ~ 1 + Group + Age + IQ + (1|SiteID)';

lme = fitlme(T, formula, 'FitMethod','REML');

%% -------------------------
% 6) 输出结果（屏幕显示 + 提取Group_ASD行）
% --------------------------
disp(lme);
anovaTbl = anova(lme, 'DFMethod','Residual');
disp(anovaTbl);

coef = lme.Coefficients;

row = strcmp(coef.Name, 'Group_ASD');
if ~any(row)
     warning('Could not find coefficient named Group_ASD. Available names:');
    disp(coef.Name);
end

%% -------------------------
% 7) 保存：模型对象 + 精简结果 + CSV + 输入快照
% --------------------------
% 7.1 保存完整模型对象（最重要）
lme_degree = lme;
save('LMM_degree_OFC.mat', 'lme_degree');

% 7.2 保存关键结果（用于文章/回复信）
Degree_LMM_Result = struct();
Degree_LMM_Result.Metric   = 'Degree centrality';
Degree_LMM_Result.ROI      = roi_OFC;
Degree_LMM_Result.Formula  = formula;
Degree_LMM_Result.N        = height(T);

if any(row)
    Degree_LMM_Result.Effect    = 'ASD vs TD (Group_ASD)';
    Degree_LMM_Result.Estimate  = coef.Estimate(row);
    Degree_LMM_Result.SE        = coef.SE(row);
    Degree_LMM_Result.tStat     = coef.tStat(row);
    Degree_LMM_Result.DF        = coef.DF(row);
    Degree_LMM_Result.pValue    = coef.pValue(row);
    Degree_LMM_Result.CI_Lower  = coef.Lower(row);
    Degree_LMM_Result.CI_Upper  = coef.Upper(row);
end

save('LMM_degree_OFC_results.mat', 'Degree_LMM_Result', 'coef', 'anovaTbl');

% 7.3 导出CSV（方便粘贴/补表）
DegreeTable = struct2table(Degree_LMM_Result);
writetable(DegreeTable, 'LMM_degree_OFC_results.csv');

% 7.4 保存输入数据
save('LMM_degree_OFC_workspace_snapshot.mat', 'T', 'y', 'Group', 'Age', 'IQ', 'SiteID', 'roi_OFC');

fprintf('\nDONE. Saved:\n');
fprintf('  LMM_degree_OFC.mat\n');
fprintf('  LMM_degree_OFC_results.mat\n');
fprintf('  LMM_degree_OFC_results.csv\n');
fprintf('  LMM_degree_OFC_workspace_snapshot.mat\n');
