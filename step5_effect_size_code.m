y_acp  = Net_para.bin.acp(:);
y_alp  = Net_para.bin.alp(:);
y_acpr = Net_para.bin.acpratio(:);
y_alpr = Net_para.bin.alpratio(:);

ES_acp  = covariate_adjusted_d(y_acp,  group, cov);
ES_alp  = covariate_adjusted_d(y_alp,  group, cov);
ES_acpr = covariate_adjusted_d(y_acpr, group, cov);
ES_alpr = covariate_adjusted_d(y_alpr, group, cov);

ES_acp
ES_alp
ES_acpr
ES_alpr



roi = 26;

y_deg   = Node_para.bin.adeg(:, roi);
y_bw    = Node_para.bin.abw(:, roi);
y_prank = Node_para.bin.aprank(:, roi);

ES_deg   = covariate_adjusted_d(y_deg,   group, cov);
ES_bw    = covariate_adjusted_d(y_bw,    group, cov);
ES_prank = covariate_adjusted_d(y_prank, group, cov);

ES_deg
ES_bw
ES_prank



edge_summary = mean(JS_26(:, ind), 2);
size(edge_summary)



ES_edge = covariate_adjusted_d(edge_summary, group, cov);
ES_edge




Metric = {
    'Global_acp'
    'Global_alp'
    'Global_acpratio'
    'Global_alpratio'
    'Nodal_ROI26_degree'
    'Nodal_ROI26_betweenness'
    'Nodal_ROI26_pagerank'
    'Edge_ROI26_65edges_summary'
    };

Cohens_d = [
    ES_acp.d
    ES_alp.d
    ES_acpr.d
    ES_alpr.d
    ES_deg.d
    ES_bw.d
    ES_prank.d
    ES_edge.d
    ];

Hedges_g = [
    ES_acp.g
    ES_alp.g
    ES_acpr.g
    ES_alpr.g
    ES_deg.g
    ES_bw.g
    ES_prank.g
    ES_edge.g
    ];

T_ES = table(Metric, Cohens_d, Hedges_g);
T_ES

writetable(T_ES,'EffectSize_CovAdjusted.csv');
save('EffectSize_CovAdjusted.mat','T_ES');
