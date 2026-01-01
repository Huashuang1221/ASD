%% network calculation
cd D:\huashuang
[Net_para, Node_para] = gretna_batch_networkanalysis_bin('D:\huashuang\ct_210_total', 'Har', 0.03, 0.4, 0.02, 100, 'pos', 's');
cd D:\huashuang
save Net_para.mat Net_para Node_para
clear

%% statistical analysis
ind_asm = 1:207;
ind_hcs = 208:424;

cd D:\huashuang
load AGE_IQ_424.mat
load Net_para.mat

[T.sw, P.sw] = gretna_permutation_ttest([Net_para.bin.acp(ind_asm) Net_para.bin.alp(ind_asm) Net_para.bin.acpratio(ind_asm) Net_para.bin.alpratio(ind_asm)], [Net_para.bin.acp(ind_hcs) Net_para.bin.alp(ind_hcs) Net_para.bin.acpratio(ind_hcs) Net_para.bin.alpratio(ind_hcs)], 10000, AGE_IQ_424(ind_asm,:), AGE_IQ_424(ind_hcs,:));

[T.acp,    P.acp]    = gretna_permutation_ttest(Node_para.bin.acp(ind_asm,:),    Node_para.bin.acp(ind_hcs,:),    10000, AGE_IQ_424(ind_asm,:), AGE_IQ_424(ind_hcs,:));
[T.adeg,   P.adeg]   = gretna_permutation_ttest(Node_para.bin.adeg(ind_asm,:),   Node_para.bin.adeg(ind_hcs,:),   10000, AGE_IQ_424(ind_asm,:), AGE_IQ_424(ind_hcs,:));
[T.age,    P.age]    = gretna_permutation_ttest(Node_para.bin.age(ind_asm,:),    Node_para.bin.age(ind_hcs,:),    10000, AGE_IQ_424(ind_asm,:), AGE_IQ_424(ind_hcs,:));
[T.abw,    P.abw]    = gretna_permutation_ttest(Node_para.bin.abw(ind_asm,:),    Node_para.bin.abw(ind_hcs,:),    10000, AGE_IQ_424(ind_asm,:), AGE_IQ_424(ind_hcs,:));
[T.aeigv,  P.aeigv]  = gretna_permutation_ttest(Node_para.bin.aeigv(ind_asm,:),  Node_para.bin.aeigv(ind_hcs,:),  10000, AGE_IQ_424(ind_asm,:), AGE_IQ_424(ind_hcs,:));
[T.aprank, P.aprank] = gretna_permutation_ttest(Node_para.bin.aprank(ind_asm,:), Node_para.bin.aprank(ind_hcs,:), 10000, AGE_IQ_424(ind_asm,:), AGE_IQ_424(ind_hcs,:));

cd D:\huashuang
save Diff_net.mat T P
clear

%% statistical analysis of ID = 26 ROI
cd D:\huashuang\ct_210_total
mats = ls('Har*.mat');

for i = 1:size(mats,1)
    load(mats(i,:))
    JS_26(i,:) = Har_data(26,:);
end
JS_26(:,26) = [];

cd D:\huashuang
load AGE_IQ_424.mat

ind_asm = 1:207;
ind_hcs = 208:424;
[T.roi26, P.roi26] = gretna_permutation_ttest(JS_26(ind_asm,:), JS_26(ind_hcs,:), 10000, AGE_IQ_424(ind_asm,:), AGE_IQ_424(ind_hcs,:));

save Diff_js_roi26.mat T P JS_26
clear

%% correlation
cd D:\huashuang
load Net_para.mat
load ADOS.mat
load AGE_IQ_424.mat

load Diff_js_roi26.mat
[Pthr,~] = gretna_FDR(P.roi26,0.05);
Ind      = P.roi26 <= Pthr;
clear P

tmp = gretna_glm(Net_para.bin.acpratio(ID),[ADOS_Data AGE_IQ_424(ID,:)],'t',1);
T.acpration = tmp.t;
P.acpration = tmp.p;
tmp = gretna_glm(Node_para.bin.adeg(ID,26),[ADOS_Data AGE_IQ_424(ID,:)],'t',1);
T.adeg26 = tmp.t;
P.adeg26 = tmp.p;
tmp = gretna_glm(Node_para.bin.abw(ID,26),[ADOS_Data AGE_IQ_424(ID,:)],'t',1);
T.abw26 = tmp.t;
P.abw26 = tmp.p;
tmp = gretna_glm(Node_para.bin.aprank(ID,26),[ADOS_Data AGE_IQ_424(ID,:)],'t',1);
T.aprank26 = tmp.t;
P.aprank26 = tmp.p;
tmp = gretna_glm(JS_26(ID,Ind),[ADOS_Data AGE_IQ_424(ID,:)],'t',1);
T.js26 = tmp.t;
P.js26 = tmp.p;

save Clinical_correlation.mat T P
clear

%% pls
cd D:\huashuang
load ADOS.mat
load AGE_IQ_424.mat

load Diff_js_roi26.mat
[Pthr,~]    = gretna_FDR(P.roi26,0.05);
Ind         = find(P.roi26 <= Pthr);
tmp         = gretna_glm(JS_26(ID,Ind), AGE_IQ_424(ID,:),'r',1);

[~,~,XS,YS,~,~,~,stats] = plsregress(tmp.r,ADOS_Data,1);
[r, p]      = corr(XS,YS);
% [r1, p1]      = corr(XS,tmp.r);





