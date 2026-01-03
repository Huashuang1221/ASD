%% =========================
% Bootstrap reproducibility (90% subsampling, 10,000 iterations)
% Output: CV for each metric
% =========================

clear; clc;
rng(1);  % for reproducibility

%% --------- User inputs ---------
B = 10000;      % bootstrap iterations
frac = 0.90;    % 90% participants each iteration

% ====== Replace these with your variables ======
% GlobalMat: Nsub x Kglobal
% NodalMat : Nsub x Knodal (optional, can be [])
GlobalMat = GlobalMat;  % <-- replace
NodalMat  = [];         % <-- replace or keep []

% Optional metric names
if exist('global_names','var') ~= 1
    global_names = arrayfun(@(i) sprintf('Global_%02d',i), 1:size(GlobalMat,2), 'UniformOutput', false);
end
if ~isempty(NodalMat) && exist('nodal_names','var') ~= 1
    nodal_names = arrayfun(@(i) sprintf('Nodal_%04d',i), 1:size(NodalMat,2), 'UniformOutput', false);
end

%% --------- Checks ---------
N = size(GlobalMat,1);
assert(N > 5, 'Too few subjects.');
n_sub = max(2, round(frac*N));

%% --------- Bootstrap: mean within subsample ---------
K1 = size(GlobalMat,2);
boot_mean_global = zeros(B, K1);

for b = 1:B
    idx = randsample(N, n_sub, false);        % 90% without replacement
    boot_mean_global(b,:) = mean(GlobalMat(idx,:), 1, 'omitnan');
end

% CV across iterations for each metric
mu_g = mean(boot_mean_global, 1, 'omitnan');
sd_g = std(boot_mean_global, 0, 1, 'omitnan');
cv_g = sd_g ./ abs(mu_g);     % use abs to avoid sign issues if mean < 0
cv_g_pct = cv_g * 100;

Tg = table(global_names(:), mu_g(:), sd_g(:), cv_g(:), cv_g_pct(:), ...
    'VariableNames', {'Metric','BootMean','BootSD','CV','CV_percent'});

writetable(Tg, 'Bootstrap_CV_Global.csv');
disp('Saved: Bootstrap_CV_Global.csv');

%% --------- Optional: nodal metrics (could be many columns) ---------
if ~isempty(NodalMat)
    assert(size(NodalMat,1) == N, 'NodalMat must have same rows as GlobalMat');
    K2 = size(NodalMat,2);
    boot_mean_nodal = zeros(B, K2);

    for b = 1:B
        idx = randsample(N, n_sub, false);
        boot_mean_nodal(b,:) = mean(NodalMat(idx,:), 1, 'omitnan');
    end

    mu_n = mean(boot_mean_nodal, 1, 'omitnan');
    sd_n = std(boot_mean_nodal, 0, 1, 'omitnan');
    cv_n = sd_n ./ abs(mu_n);
    cv_n_pct = cv_n * 100;

    Tn = table(nodal_names(:), mu_n(:), sd_n(:), cv_n(:), cv_n_pct(:), ...
        'VariableNames', {'Metric','BootMean','BootSD','CV','CV_percent'});

    writetable(Tn, 'Bootstrap_CV_Nodal.csv');
    disp('Saved: Bootstrap_CV_Nodal.csv');
end
