r_mol  = results.mean_r(:);         % 26x1 effect size (spatial r)
p_spin = results.mean_spin_p(:);    % 26x1 spin-test p

% map names
if isfield(results,'names')
    mol_names = results.names(:);
else
    mol_names = arrayfun(@(i) sprintf('Map_%02d', i), (1:numel(r_mol))', 'UniformOutput', false);
end

% ===== 1) within-family FDR (BH) on spin p =====
q_fdr = mafdr(p_spin,'BHFDR',true);
sig_fdr = q_fdr < 0.05;

fprintf('Molecular family: n_maps=%d, n_sig(q<0.05)=%d\n', numel(p_spin), sum(sig_fdr));

% ===== 2) Bootstrap Z (for Fig 4E–F) =====

if isfield(results,'trans_weight_boot_z')
    z_boot = results.trans_weight_boot_z(:);
elseif isfield(results,'z_boot')
    z_boot = results.z_boot(:);
else
    warning('No bootstrap Z field found. Please set z_boot = your bootstrap-derived Z-scores (26x1).');
    z_boot = nan(size(r_mol));
end

% ===== 3) Export supplement table =====
MolTable = table( ...
    mol_names, r_mol, r_mol.^2, p_spin, q_fdr, z_boot, sig_fdr, ...
    'VariableNames', {'NeuroMap','r','r2','p_spin','q_fdr','z_boot','sig_fdr'} );

writetable(MolTable, 'SupTable_Molecular_26maps.csv');

% ===== 4) Define family-level primary p for across-family correction =====

if isfield(results,'p_spin_pls1')
    p_primary_mol = results.p_spin_pls1;
    primary_def = 'PLS1 multivariate correspondence spin-test p (model-level)';
else
    p_primary_mol = min(p_spin); % fallback
    primary_def = 'minimum spin-test p across 26 maps (fallback; not model-level)';
end

% ===== 5) Save MolecularResults for Bonferroni_Holm =====
MolecularResults = struct();
MolecularResults.family = 'Molecular';
MolecularResults.map_names = mol_names;
MolecularResults.r = r_mol;
MolecularResults.p_spin = p_spin;
MolecularResults.q_fdr = q_fdr;
MolecularResults.z_boot = z_boot;
MolecularResults.sig_fdr = sig_fdr;

MolecularResults.p_primary = p_primary_mol;
MolecularResults.primary_definition = primary_def;

save('Results_Molecular.mat','MolecularResults');

disp('Done. Saved: Results_Molecular.mat and SupTable_Molecular_26maps.csv');