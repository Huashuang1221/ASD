%% corr the trans and z
%% para
rng(0)
load('Diff_js_roi26.mat')
data = T.roi26.real';
data = [data(1:25);nan;data(26:end)];
data_all = nan(210,1);
data_all(setdiff(1:210,[59 163])) = data;
data = data_all;
load('spin_per_Fan210.mat')

load('neuro_trans_net_Fan_210.mat')
TM_data = TM_data.Fan_210;

pls_dim = 1;
spin_time = 10000;
boot_time = 10000;
p_thr = 0.05;

results = struct;


%% cal
% corr
[results.mean_r,results.mean_p] = corr(TM_data,data,'rows','pairwise','type','Spearman');
mean_r_spin = zeros(size(TM_data,2),spin_time);
for i_spin = 1:spin_time
    mean_r_spin(:,i_spin) = corr(TM_data,data(perm_id(:,i_spin)),'rows','pairwise','type','Spearman');
end
results.mean_spin_p = (sum(abs(mean_r_spin) >= repmat(abs(results.mean_r),1,spin_time),2) + 1)/(spin_time + 1);
results.mean_spin_fdr = gretna_FDR(results.mean_spin_p,p_thr);
if ~isempty(results.mean_spin_fdr)
    sig_index = find(results.mean_spin_p <= results.mean_spin_fdr);
    results.mean_spin_sig = [trans_name(sig_index),...
        num2cell(results.mean_r(sig_index)),...
        num2cell(results.mean_spin_p(sig_index))];
end

% pls
del_region = isnan(sum(TM_data,2)) | isnan(data);
[~,~,XS,YS,~,PCTVAR,~,stats] = plsregress(TM_data(~del_region,:),data(~del_region),pls_dim);

results.YS = YS;
results.PCTVAR = PCTVAR;
results.X_Score = XS;
results.trans_weight = stats.W;

% spin
results.PCTVAR_spin = zeros(spin_time,pls_dim);
for i_spin = 1:spin_time
    data_spin = data(perm_id(:,i_spin));
    del_region = isnan(sum(TM_data,2)) | isnan(data_spin);
    [~,~,~,~,~,PCTVAR,~,~] = plsregress(TM_data(~del_region,:),data_spin(~del_region),pls_dim);
    results.PCTVAR_spin(i_spin,:) = PCTVAR(2,:);
end
results.p_spin = (sum(results.PCTVAR_spin >= repmat(results.PCTVAR(2,:),spin_time,1))+1)/(spin_time + 1);

% boot
W_boot = zeros(length(results.trans_weight),boot_time);
for iboot = 1:boot_time
    if mod(iboot,boot_time/10) == 0
        disp([num2str(iboot) '/' num2str(boot_time) ')  ' datestr(clock)]);
    end

    myresample = randsample(length(data),length(data),1);
    X_boot = TM_data(myresample,:);
    Y_boot = data(myresample);

    del_region = isnan(sum(X_boot,2)) | isnan(Y_boot);
    [~,~,~,~,~,~,~,stats] = plsregress(X_boot(~del_region,:),Y_boot(~del_region),pls_dim);
    W_boot(:,iboot) = stats.W;
end
results.trans_weight_boot = W_boot;
W_boot = std(W_boot,[],2);

results.trans_weight_boot_z = results.trans_weight./W_boot;
results.trans_weight_boot_p = 2 * (1 - normcdf(abs(results.trans_weight_boot_z)));

results.trans_sig = [trans_name(results.trans_weight_boot_p <= p_thr),...
    num2cell(results.trans_weight(results.trans_weight_boot_p <= p_thr)),...
    num2cell(results.trans_weight_boot_z(results.trans_weight_boot_p <= p_thr)),...
    num2cell(results.trans_weight_boot_p(results.trans_weight_boot_p <= p_thr))];

save('results_corr_trans.mat','results')