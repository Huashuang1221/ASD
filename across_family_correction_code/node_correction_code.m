size(Node_para.bin.adeg)
size(Node_para.bin.abw)
size(Node_para.bin.aprank)

roi_OFC = 26;
idx_ASD = 1:207;
idx_TD  = 208:424;

metrics_auc = {'adeg','abw','aprank','aeigv','age'}; % 你要哪些就留哪些（字段必须存在）

for i = 1:numel(metrics_auc)
    met = metrics_auc{i};
    X = Node_para.bin.(met);     % 424×208
    v = X(:, roi_OFC);           % 424×1

    x_ASD = v(idx_ASD);
    x_TD  = v(idx_TD);

    [~, p_raw, ~, stats] = ttest2(x_ASD, x_TD);

    sp = sqrt(((numel(x_ASD)-1)*var(x_ASD) + (numel(x_TD)-1)*var(x_TD)) / (numel(x_ASD)+numel(x_TD)-2));
    d  = (mean(x_ASD) - mean(x_TD)) / sp;

    fprintf('%-8s ROI%03d: t=%.2f, df=%d, p=%.6g, d=%.2f\n', met, roi_OFC, stats.tstat, stats.df, p_raw, d);
end


% 把5个 raw p 收集起来（顺序你自己定，但要和指标名一致）
p_nodal = [ ...
    7.38823e-05;   % adeg
    5.05053e-05;   % abw
    8.53353e-05;   % aprank
    4.60081e-04;   % aeigv
    1.17307e-04];  % age

% BH-FDR 校正（family 内）
q_nodal = mafdr(p_nodal, 'BHFDR', true);

% 看哪些 survives
[q_nodal, q_nodal < 0.05]


format long g
[q_nodal, q_nodal<0.05]


format long g

metric = {'adeg','abw','aprank','aeigv','age'}';

tval = [-4.00; -4.10; -3.97; -3.53; -3.89];
df   = [422; 422; 422; 422; 422];

p_raw = p_nodal(:);
q_fdr = q_nodal(:);

dval = [-0.39; -0.40; -0.39; -0.34; -0.38];  % 你按打印结果填；或把循环里 d 存起来

T = table(metric, tval, df, p_raw, q_fdr, dval);

disp(T)

% 导出补充表
writetable(T, 'SupTable_OFC_Nodal_FDR.csv');


p_primary_nodal = 7.38823e-05;  % adeg 的 raw p

NodalResults = struct();

NodalResults.family = 'Nodal_OFC';

NodalResults.metric = metric;        % {'adeg','abw',...}
NodalResults.tval   = tval;
NodalResults.df     = df;
NodalResults.p_raw  = p_raw;
NodalResults.q_fdr  = q_fdr;
NodalResults.dval   = dval;

% —— 最重要的两行 —— %
NodalResults.primary_metric = 'adeg';
NodalResults.p_primary      = p_primary_nodal;  % 7.38823e-05

save('Results_Nodal_OFC.mat', 'NodalResults');


