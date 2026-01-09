cd D:\huashuang
load('bin_65_Diff_js_roi26.mat');   % JS_26, P (and maybe T, but we won't use T)

% ===== 0) 基本信息 =====
idx_ASD = 1:207;
idx_TD  = 208:424;
nROI = 208;

% ===== 1) edge raw p (来自置换检验，已控制协变量) =====
p_edge = P.roi26(:);                % 207×1
q_edge = mafdr(p_edge,'BHFDR',true);
sig    = q_edge < 0.05;

fprintf('Edge family: n_edges=%d, n_sig(q<0.05)=%d\n', numel(p_edge), sum(sig));

% ===== 2) 映射到真实 target ROI 编号 =====
roi_list   = setdiff(1:nROI, 26);
target_roi = roi_list(:);           % 207×1

% ===== 3) 重算每条边的 t 值（用于描述；p值仍用置换p） =====
t_edge = nan(207,1);
for e = 1:207
    [~,~,~,stat] = ttest2(JS_26(idx_ASD,e), JS_26(idx_TD,e), 'Vartype','unequal');
    t_edge(e) = stat.tstat;
end

% ===== 4) Cohen's d（每条边） =====
d_edge = nan(207,1);
for e = 1:207
    xA = JS_26(idx_ASD, e);
    xT = JS_26(idx_TD,  e);
    sp = sqrt(((numel(xA)-1)*var(xA) + (numel(xT)-1)*var(xT)) / (numel(xA)+numel(xT)-2));
    d_edge(e) = (mean(xA) - mean(xT)) / sp;
end

% ===== 5) 导出补充表 =====
EdgeTable_all = table( ...
    repmat(26,207,1), target_roi, t_edge, p_edge, q_edge, d_edge, sig, ...
    'VariableNames', {'roi_seed','roi_target','tval','p_raw','q_fdr','cohen_d','sig_fdr'});

writetable(EdgeTable_all, 'SupTable_Edge_ROI26_all207.csv');

EdgeTable_sig = EdgeTable_all(sig,:);
writetable(EdgeTable_sig, 'SupTable_Edge_ROI26_sig.csv');

% ===== 6) 为 across-family 准备 edge primary p =====
p_primary_edge = min(p_edge(sig));  % 推荐：在FDR显著边中取最小raw p

% ===== 7) 保存 EdgeResults（后续跨-family一键调用） =====
EdgeResults = struct();
EdgeResults.family = 'Edge_OFC_ROI26';
EdgeResults.roi_seed = 26;
EdgeResults.roi_target = target_roi;
EdgeResults.tval = t_edge;          % 描述性t
EdgeResults.p_raw = p_edge;         % 置换p（控制协变量）
EdgeResults.q_fdr = q_edge;
EdgeResults.cohen_d = d_edge;
EdgeResults.sig_fdr = sig;

EdgeResults.p_primary = p_primary_edge;
EdgeResults.primary_definition = 'min raw p among FDR-significant edges (q<0.05)';
EdgeResults.note = 'p_raw from permutation test with covariates; tval computed by Welch t-test for descriptive reporting';

save('Results_Edge_OFC.mat','EdgeResults');
