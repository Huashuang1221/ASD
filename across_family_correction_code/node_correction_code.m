
% Prerequisites:
% 1) load Node_para.mat
% 2) load AGE_IQ_424.mat 
% 3) ASD = 1:207, TD = 208:424  

roi_OFC = 26;
idx_ASD = 1:207;
idx_TD  = 208:424;

% ---- covariates (age + IQ) ----
cov = AGE_IQ_424(:, [1 2]);  % adjust columns if needed

% ---- metrics to test (must exist in Node_para.bin) ----
metrics = {'adeg','abw','aprank','aeigv','age'}; 


thr = 0.03:0.02:0.4;  

nM = numel(metrics);
tval  = nan(nM,1);
df    = nan(nM,1);
p_raw = nan(nM,1);
dval  = nan(nM,1);

for i = 1:nM
    met = metrics{i};
    X = Node_para.bin.(met);

    % ---- extract subject-level AUC value for OFC node ----
    if ndims(X) == 2
        % assume: sub × node
        v = X(:, roi_OFC);
    elseif ndims(X) == 3
        % assume: thr × sub × node  (common in GRETNA outputs)
        v = squeeze(trapz(thr, X(:,:,roi_OFC), 1));  % 1×sub
        v = v(:);
    else
        error('Unexpected dimension for Node_para.bin.%s', met);
    end

    % ---- residualize (remove age + IQ) ----
    v_res = residualize_covariates(v, cov);

    % ---- group comparison on residuals ----
    [~, p, ~, stats] = ttest2(v_res(idx_ASD), v_res(idx_TD));
    p_raw(i) = p;
    tval(i)  = stats.tstat;
    df(i)    = stats.df;

    % ---- effect size on residuals ----
    dval(i) = cohens_d(v_res(idx_ASD), v_res(idx_TD));
end

% ---- within-family BH-FDR ----
q_fdr = mafdr(p_raw, 'BHFDR', true);

% ---- export Supplement table ----
T = table(metrics(:), tval, df, p_raw, q_fdr, dval, ...
    'VariableNames', {'metric','t','df','p_raw','q_fdr','cohens_d'});
writetable(T, 'SupTable_OFC_Nodal_FDR.csv');

% ---- cross-family Holm: define a pre-specified primary metric ----
primary_metric = 'adeg'; % e.g., degree centrality (predefined)
primary_idx = find(strcmp(metrics, primary_metric), 1);
p_primary_nodal = p_raw(primary_idx);

% ---- save for downstream Holm script ----
NodalResults = struct();
NodalResults.family = 'Nodal_OFC';
NodalResults.primary_metric = primary_metric;
NodalResults.p_primary = p_primary_nodal;
NodalResults.table = T;
save('Results_Nodal_OFC.mat', 'NodalResults');

disp(T);

%% -------- local helper functions (MATLAB R2018b supports local functions in scripts) --------
function resid = residualize_covariates(y, cov)
    y = y(:);
    X = [ones(size(cov,1),1), cov];
    beta = X \ y;
    resid = y - X * beta;
end

function d = cohens_d(x1, x2)
    x1 = x1(:); x2 = x2(:);
    n1 = numel(x1); n2 = numel(x2);
    s1 = var(x1, 1);
    s2 = var(x2, 1);
    sp = sqrt(((n1-1)*s1 + (n2-1)*s2) / (n1 + n2 - 2));
    d = (mean(x1) - mean(x2)) / sp;
end

