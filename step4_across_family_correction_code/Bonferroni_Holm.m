%% Bonferroni_Holm_UPDATED.m
% Across-family Bonferroni–Holm step-down correction on family-level primary p-values

load('Results_Global.mat');      % GlobalResults
load('Results_Nodal_OFC.mat');   % NodalResults
load('Results_Molecular.mat');   % MolecularResults
load('Results_Edge_OFC.mat');    % EdgeResults

% --- Get family-level primary p (one per family) ---
% Global：
if isfield(GlobalResults,'p_primary')
    p_global = GlobalResults.p_primary;
    def_global = GlobalResults.primary_definition;
else
    % fallback: use raw p (but please standardize later)
    p_global = GlobalResults.p_raw;
    def_global = 'fallback: GlobalResults.p_raw';
end

p_nodal = NodalResults.p_primary;

% Edge
p_edge  = EdgeResults.p_primary;

% neurotransmitter
p_mol   = MolecularResults.p_primary;

p_raw = [p_global, p_nodal, p_edge, p_mol];
labels = {'Global','Nodal_OFC','Edge_OFC','Molecular'}';

%% Holm step-down
m = numel(p_raw);
[p_sorted, idx] = sort(p_raw,'ascend');
adj_sorted = zeros(size(p_sorted));

for k = 1:m
    adj_sorted(k) = (m-k+1) * p_sorted(k);
end

adj_sorted = cummax(adj_sorted);
adj_sorted(adj_sorted>1) = 1;

p_holm = zeros(size(p_raw));
p_holm(idx) = adj_sorted;

%% Display + Save
T = table(labels, p_raw(:), p_holm(:), 'VariableNames', {'Family','p_primary','p_Holm'});
disp(T)

save('AcrossFamily_BonferroniHolm.mat','T')
writetable(T,'Table_S_AcrossFamily_BonferroniHolm.csv');
