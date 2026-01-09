%% =========================
% Weighted MBN (robustness check)
% =========================

Mat_path = 'D:\huashuang\ct_210_total';  
File_filter = 'Har';                    

Thr1 = 0.03; Thr2 = 0.40; Delta = 0.02;
N_rand = 100;
Stype = 'pos';
Ttype = 's';


[Net_para_w, Node_para_w] = gretna_batch_networkanalysis(Mat_path, File_filter, Thr1, Thr2, Delta, N_rand, Stype, Ttype);

save Net_para_weighted.mat Net_para_w Node_para_w
