load OFC_adeg.mat            % 变量名假设叫 OFC_adeg（如果不是，whos 看一下）
load MEDICINE_ID.mat         % 变量名假设叫 MEDICINE_ID
load non_MEDICINE_ID.mat     % 变量名假设叫 non_MEDICINE_ID

% 1) 取 ASD 组的 OFC degree（第26列）
ofc_asd = OFC_adeg(1:207, 26);     % 207×1

% 2) 用“序号索引”直接提取
ofc_medicated = ofc_asd(MEDICINE_ID);        % 用药 ASD 的 OFC degree
ofc_drugnaive = ofc_asd(non_MEDICINE_ID);    % 未用药 ASD 的 OFC degree

% 3) 取 TD 的 OFC degree（之后算效应量要用）
ofc_td = OFC_adeg(208:424, 26);    % 217×1

% 确保索引都在 1~207
assert(all(MEDICINE_ID>=1 & MEDICINE_ID<=207), 'MEDICINE_ID 超出 1~207');
assert(all(non_MEDICINE_ID>=1 & non_MEDICINE_ID<=207), 'non_MEDICINE_ID 超出 1~207');

% 检查两组是否有重叠
overlap = intersect(MEDICINE_ID, non_MEDICINE_ID);
assert(isempty(overlap), '用药与未用药ID有重叠，请检查！');

% 检查是否覆盖全部ASD（可选）
n_missing = 207 - numel(unique([MEDICINE_ID(:); non_MEDICINE_ID(:)]));
fprintf('ASD中用药=%d, 未用药=%d, 未分类/缺失=%d\n', numel(MEDICINE_ID), numel(non_MEDICINE_ID), n_missing);

MEDICINE_ID     = unique(MEDICINE_ID(:));
non_MEDICINE_ID = unique(non_MEDICINE_ID(:));



% 计算 Cohen's d 的匿名函数（两独立样本）
cohens_d = @(x,y) (mean(x,'omitnan') - mean(y,'omitnan')) ./ ...
    sqrt(((numel(x)-1)*var(x,'omitnan') + (numel(y)-1)*var(y,'omitnan')) / (numel(x)+numel(y)-2));

d_naive = cohens_d(ofc_drugnaive, ofc_td);
d_medi  = cohens_d(ofc_medicated, ofc_td);

fprintf('drug-naive ASD vs TD: n=%d vs %d, d=%.3f\n', numel(ofc_drugnaive), numel(ofc_td), d_naive);
fprintf('medicated ASD vs TD: n=%d vs %d, d=%.3f\n', numel(ofc_medicated), numel(ofc_td), d_medi);
