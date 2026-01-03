cd D:\AAARESULTS_ASD_SMRI\wang_jinhui_results\jinhui_rsults_code
gamma_mat = Net_para.bin.cpratio;   % 19 × 424

% ⚠️ 把这里改成你真实用的阈值范围
Thr1 = 0.10;
Thr2 = 0.28;
Delta = 0.01;

thr = Thr1:Delta:Thr2;

% 简单校验
assert(numel(thr) == size(gamma_mat,1), '阈值数量与19不一致，请检查 Thr1/Thr2/Delta');

auc_gamma = trapz(thr, gamma_mat) / (Thr2 - Thr1);  % 1×424
auc_gamma = auc_gamma(:);                            % 424×1


idx_ASD = 1:207;
idx_TD  = 208:424;

x_ASD = auc_gamma(idx_ASD);
x_TD  = auc_gamma(idx_TD);

[~, p_raw, ~, stats] = ttest2(x_ASD, x_TD);

t_raw  = stats.tstat;
df_raw = stats.df;

% Cohen's d
sp = sqrt(((numel(x_ASD)-1)*var(x_ASD) + (numel(x_TD)-1)*var(x_TD)) / ...
          (numel(x_ASD)+numel(x_TD)-2));
d_raw = (mean(x_ASD) - mean(x_TD)) / sp;

fprintf('Normalized Cp (gamma, AUC): t = %.2f, df = %d, p = %.6g, d = %.2f\n', ...
        t_raw, df_raw, p_raw, d_raw);
