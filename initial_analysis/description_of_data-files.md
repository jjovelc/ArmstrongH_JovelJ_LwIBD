# Description of Data Files

This document provides a comprehensive description of all files in the project directory structure.

<table>
<thead>
<tr>
<th>Directory</th>
<th>File Name</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr>
<td rowspan="4">Root Directory</td>
<td>all_samples_metaph3-estCounts.tsv</td>
<td>Metaphlan3 estimated counts for all samples</td>
</tr>
<tr>
<td>all_samples_metaph3-estCounts_shortNames.tsv</td>
<td>Metaphlan3 estimated counts with shortened names for all samples</td>
</tr>
<tr>
<td>kraken2_all_samples_mpa_tax.tsv</td>
<td>Kraken2 taxonomic classification results for all samples in MPA format</td>
</tr>
<tr>
<td>kraken2_all_samples_mpa_tax_shortNames.tsv</td>
<td>Kraken2 taxonomic classification with shortened names for all samples</td>
</tr>
<tr>
<td rowspan="4">Root Directory</td>
<td>metadata.tsv</td>
<td>Metadata file containing sample information (sex, severity, diagnosis, timePoint, patient)</td>
</tr>
<tr>
<td>metabolic_profiling_rxn_estAbundance.tsv</td>
<td>Metabolic profiling reaction estimated abundance data</td>
</tr>
<tr>
<td>flare-wk0_vs_remission-wk52_pAdj0.05.tsv</td>
<td>Differential analysis results comparing flare week 0 vs remission week 52 (adjusted p-value < 0.05)</td>
</tr>
<tr>
<td>remission-wk0_vs_flare-wk52_pAdj0.05.tsv</td>
<td>Differential analysis results comparing remission week 0 vs flare week 52 (adjusted p-value < 0.05)</td>
</tr>
<tr>
<td rowspan="7">exploratory_plots</td>
<td>HC_plot.pdf</td>
<td>Hierarchical clustering plot in PDF format</td>
</tr>
<tr>
<td>PCA_plot.pdf</td>
<td>Principal Component Analysis plot in PDF format</td>
</tr>
<tr>
<td>PCoA_by_diagnosis.png</td>
<td>Principal Coordinates Analysis plot colored by diagnosis</td>
</tr>
<tr>
<td>PCoA_by_severity.png</td>
<td>Principal Coordinates Analysis plot colored by disease severity</td>
</tr>
<tr>
<td>PCoA_by_sex.png</td>
<td>Principal Coordinates Analysis plot colored by sex</td>
</tr>
<tr>
<td>PCoA_by_timePoint.png</td>
<td>Principal Coordinates Analysis plot colored by time point</td>
</tr>
<tr>
<td>PCoA_plot.pdf</td>
<td>Principal Coordinates Analysis plot in PDF format</td>
</tr>
<tr>
<td rowspan="60">0remission_52flare_analyses_230913</td>
<td>analysis_circusplot.Rmd</td>
<td>R Markdown script for creating circus plots</td>
</tr>
<tr>
<td>correlation_circosPlots.R</td>
<td>R script for generating correlation circos plots</td>
</tr>
<tr>
<td>desc_columns_allData.tsv</td>
<td>Description of columns in all data files</td>
</tr>
<tr>
<td>DESeq2_microbiome.R</td>
<td>R script for differential expression analysis using DESeq2 on microbiome data</td>
</tr>
<tr>
<td>diff_data.txt</td>
<td>Differential analysis data in text format</td>
</tr>
<tr>
<td>diff_table.txt</td>
<td>Differential analysis results table</td>
</tr>
<tr>
<td>energy_and_consumption_0-52w.tsv</td>
<td>Energy and consumption data for weeks 0-52</td>
</tr>
<tr>
<td>extract_difference_52-0.py</td>
<td>Python script to extract differences between week 52 and week 0</td>
</tr>
<tr>
<td>helper_functions.R</td>
<td>R script containing helper functions for data analysis</td>
</tr>
<tr>
<td>hmdb_metabolites.xml</td>
<td>HMDB (Human Metabolome Database) metabolites data in XML format</td>
</tr>
<tr>
<td>hmdb_metabolites.zip</td>
<td>Compressed archive of HMDB metabolites data</td>
</tr>
<tr>
<td>list_mets_metabolitesNames.txt</td>
<td>List of metabolite names</td>
</tr>
<tr>
<td>meta.tsv</td>
<td>Metadata file in TSV format</td>
</tr>
<tr>
<td>metabolites_0remission_52flare_4metaboAnalyst.txt</td>
<td>Metabolites data formatted for MetaboAnalyst (remission week 0 vs flare week 52)</td>
</tr>
<tr>
<td>metabolites_0remission_52flare_annotation-b.txt</td>
<td>Metabolite annotations (variant b) for remission week 0 vs flare week 52</td>
</tr>
<tr>
<td>metabolites_0remission_52flare_annotation.txt</td>
<td>Metabolite annotations for remission week 0 vs flare week 52</td>
</tr>
<tr>
<td>metabolites_0remission_52flare-nonNorm.txt</td>
<td>Non-normalized metabolites data for remission week 0 vs flare week 52</td>
</tr>
<tr>
<td>metabolites.txt</td>
<td>Metabolites data in text format</td>
</tr>
<tr>
<td>metabolomics_only_0-52w-b.tsv</td>
<td>Metabolomics data for weeks 0-52 (variant b, transposed)</td>
</tr>
<tr>
<td>metabolomics_only_0-52w-t.tsv</td>
<td>Metabolomics data for weeks 0-52 (transposed)</td>
</tr>
<tr>
<td>metabolomics_only_0-52w.tsv</td>
<td>Metabolomics data for weeks 0-52</td>
</tr>
<tr>
<td>metabolomics.R</td>
<td>R script for metabolomics data analysis</td>
</tr>
<tr>
<td>metadata_0-52w.tsv</td>
<td>Metadata for weeks 0-52</td>
</tr>
<tr>
<td>metadata_previous_0-52w.tsv</td>
<td>Previous version of metadata for weeks 0-52</td>
</tr>
<tr>
<td>metanb</td>
<td>Metadata notebook file</td>
</tr>
<tr>
<td>multiOmics_data_0-52w.tsv</td>
<td>Multi-omics integrated data for weeks 0-52</td>
</tr>
<tr>
<td>normalize_kraken_counts.py</td>
<td>Python script to normalize Kraken count data</td>
</tr>
<tr>
<td>parse_hmdb_xml.py</td>
<td>Python script to parse HMDB XML files</td>
</tr>
<tr>
<td>parse_xml.py</td>
<td>Python script for parsing XML files</td>
</tr>
<tr>
<td>pbmcsca_seurat_object.rds</td>
<td>Seurat object containing PBMC single-cell analysis data</td>
</tr>
<tr>
<td>results_Utest.txt</td>
<td>Results from U-test (Mann-Whitney test) analysis</td>
</tr>
<tr>
<td>rJava_1.0-11.tar.gz</td>
<td>Compressed R package archive for rJava</td>
</tr>
<tr>
<td>rowdata.txt</td>
<td>Row data information</td>
</tr>
<tr>
<td>sample_and_progression.txt</td>
<td>Sample information and disease progression data</td>
</tr>
<tr>
<td>taxa_only_0-52w_100up_CPM.tsv</td>
<td>Taxonomic data for weeks 0-52, filtered to taxa with 100+ counts, in CPM (counts per million) format</td>
</tr>
<tr>
<td>taxa_only_0-52w_100up.tsv</td>
<td>Taxonomic data for weeks 0-52, filtered to taxa with 100+ counts</td>
</tr>
<tr>
<td>taxa_only_0-52w_RF_counts_allRes_normData.tsv</td>
<td>All DESeq2 results with normalized data for RF (remission/flare) counts</td>
</tr>
<tr>
<td>taxa_only_0-52w_RF_counts_HCheatmap.pdf</td>
<td>Hierarchical clustering heatmap PDF for RF counts</td>
</tr>
<tr>
<td>taxa_only_0-52w_RF_counts_p0.01_DE_features.tsv</td>
<td>Differentially expressed features with p-value < 0.01 for RF counts</td>
</tr>
<tr>
<td>taxa_only_0-52w_RF_counts_p0.05_DE_features.tsv</td>
<td>Differentially expressed features with p-value < 0.05 for RF counts</td>
</tr>
<tr>
<td>taxa_only_0-52w_RF_counts_p0.1_DE_features.tsv</td>
<td>Differentially expressed features with p-value < 0.1 for RF counts</td>
</tr>
<tr>
<td>taxa_only_0-52w_RF_counts_PCoA.pdf</td>
<td>Principal Coordinates Analysis plot PDF for RF counts</td>
</tr>
<tr>
<td>taxa_only_0-52w_RF_counts_sign_DESeq2_res.xlsx</td>
<td>Significant DESeq2 results in Excel format for RF counts</td>
</tr>
<tr>
<td>taxa_only_0-52w_RF_counts_volcanoPlot.pdf</td>
<td>Volcano plot PDF for RF counts differential analysis</td>
</tr>
<tr>
<td>taxa_only_0-52w_RF_counts.tsv</td>
<td>Taxonomic count data for weeks 0-52 comparing remission and flare</td>
</tr>
<tr>
<td>taxa_only_0-52w_RF_meta.tsv</td>
<td>Metadata for RF (remission/flare) comparison</td>
</tr>
<tr>
<td>taxa_only_0-52w-t.tsv</td>
<td>Taxonomic data for weeks 0-52 (transposed)</td>
</tr>
<tr>
<td>taxa_only_0-52w.tsv</td>
<td>Taxonomic data for weeks 0-52</td>
</tr>
<tr>
<td>Ttest_results_metabolomics_remission0-flare52.tsv</td>
<td>T-test results for metabolomics comparing remission week 0 vs flare week 52</td>
</tr>
<tr>
<td>Untitled.numbers</td>
<td>Apple Numbers spreadsheet file (untitled)</td>
</tr>
<tr>
<td>Utest_RR_RF_0-52.R</td>
<td>R script for U-test comparing RR (remission-remission) and RF (remission-flare) for weeks 0-52</td>
</tr>
<tr>
<td>xlsx_0.6.5.tar.gz</td>
<td>Compressed R package archive for xlsx</td>
</tr>
<tr>
<td rowspan="4">0remission_52flare_analyses_230913/ewandson_example</td>
<td>analysis_circusplot.Rmd</td>
<td>R Markdown script for creating circus plots (example from Ewandson)</td>
</tr>
<tr>
<td>helper_functions.R</td>
<td>Helper functions R script (example from Ewandson)</td>
</tr>
<tr>
<td>SGCCA_HFHS.csv</td>
<td>SGCCA (Sparse Generalized Canonical Correlation Analysis) results for HFHS dataset</td>
</tr>
<tr>
<td>SGCCA_SD.csv</td>
<td>SGCCA results for SD dataset</td>
</tr>
<tr>
<td rowspan="12">0remission_52flare_analyses_230913/MetaboAnalyst_plots</td>
<td>Acetoacetate_boxplot.png</td>
<td>Boxplot visualization for Acetoacetate metabolite</td>
</tr>
<tr>
<td>acetyl_methionine_boxplots.png</td>
<td>Boxplot visualizations for acetyl methionine metabolites</td>
</tr>
<tr>
<td>correlationPlot_remission0-flare52.png</td>
<td>Correlation plot comparing remission week 0 and flare week 52</td>
</tr>
<tr>
<td>heatmap_OneFactor_MetaboAnalyst.png</td>
<td>One-factor heatmap generated by MetaboAnalyst</td>
</tr>
<tr>
<td>hippurate_boxplot.png</td>
<td>Boxplot visualization for Hippurate metabolite</td>
</tr>
<tr>
<td>metabo_VP.png</td>
<td>Volcano plot for metabolomics data</td>
</tr>
<tr>
<td>metabolites_heatmap_metaboAnalyst.png</td>
<td>Heatmap of metabolites generated by MetaboAnalyst</td>
</tr>
<tr>
<td>methionine_boxplot.png</td>
<td>Boxplot visualization for Methionine metabolite</td>
</tr>
<tr>
<td>PCA_metaboAnalyist.png</td>
<td>Principal Component Analysis plot from MetaboAnalyst</td>
</tr>
<tr>
<td>PCA_remission0_flare52.png</td>
<td>PCA plot comparing remission week 0 and flare week 52</td>
</tr>
<tr>
<td>SAM_MetaboAnalyst.png</td>
<td>SAM (Significance Analysis of Microarrays) plot from MetaboAnalyst</td>
</tr>
<tr>
<td>volcanoPlot_MetaboAnalyst.png</td>
<td>Volcano plot generated by MetaboAnalyst</td>
</tr>
<tr>
<td rowspan="64">0remission_52flare_analyses_230913/mixOmics</td>
<td>1_1_1_1_cor_matrix_circos.csv</td>
<td>Correlation matrix in CSV format for circos plot generation</td>
</tr>
<tr>
<td>ACETOACETATE_Ttest_MetaboAnalyst.png</td>
<td>T-test results visualization for Acetoacetate from MetaboAnalyst</td>
</tr>
<tr>
<td>analysis_circusplot.Rmd</td>
<td>R Markdown script for creating circus plots in mixOmics analysis</td>
</tr>
<tr>
<td>correlation_circosPlots.R</td>
<td>R script for generating correlation circos plots</td>
</tr>
<tr>
<td>correlation_matrix.tsv</td>
<td>Correlation matrix in TSV format</td>
</tr>
<tr>
<td>data</td>
<td>Data file (binary or text format)</td>
</tr>
<tr>
<td>edited_metabolites.tsv</td>
<td>Edited metabolites data in TSV format</td>
</tr>
<tr>
<td>energy_and_consumption_0-52w-b.tsv</td>
<td>Energy and consumption data for weeks 0-52 (variant b, transposed)</td>
</tr>
<tr>
<td>energy_and_consumption_0-52w.tsv</td>
<td>Energy and consumption data for weeks 0-52</td>
</tr>
<tr>
<td>found_names</td>
<td>File containing found metabolite names</td>
</tr>
<tr>
<td>HC_remission0_flare52.png</td>
<td>Hierarchical clustering plot for remission week 0 vs flare week 52</td>
</tr>
<tr>
<td>header</td>
<td>Header information file</td>
</tr>
<tr>
<td>helper_functions.R</td>
<td>Helper functions R script for mixOmics analysis</td>
</tr>
<tr>
<td>ids.selCol.tsv</td>
<td>Selected column IDs in TSV format</td>
</tr>
<tr>
<td>liverToxicity.R</td>
<td>R script for liver toxicity analysis (example dataset)</td>
</tr>
<tr>
<td>metabolites_0remission_52flare_4metaboAnalyst.txt</td>
<td>Metabolites data formatted for MetaboAnalyst</td>
</tr>
<tr>
<td>metabolites_0remission_52flare_annotation.txt</td>
<td>Metabolite annotations for remission week 0 vs flare week 52</td>
</tr>
<tr>
<td>metabolites_0remission_52flare-nonNorm.txt</td>
<td>Non-normalized metabolites data</td>
</tr>
<tr>
<td>metabolites_0remission_52flare.txt</td>
<td>Metabolites data for remission week 0 vs flare week 52</td>
</tr>
<tr>
<td>metabolites.tsv</td>
<td>Metabolites data in TSV format</td>
</tr>
<tr>
<td>metabolites.txt</td>
<td>Metabolites data in text format</td>
</tr>
<tr>
<td>metabolomics_only_0-52w_norm_remission.tsv</td>
<td>Normalized metabolomics data for weeks 0-52 (remission samples only)</td>
</tr>
<tr>
<td>metabolomics_only_0-52w_norm.tsv</td>
<td>Normalized metabolomics data for weeks 0-52</td>
</tr>
<tr>
<td>metabolomics_only_0-52w.tsv</td>
<td>Metabolomics data for weeks 0-52</td>
</tr>
<tr>
<td>metabolomics.R</td>
<td>R script for metabolomics analysis</td>
</tr>
<tr>
<td>metadata_0-52w-b.tsv</td>
<td>Metadata for weeks 0-52 (variant b, transposed)</td>
</tr>
<tr>
<td>metadata_0-52w.tsv</td>
<td>Metadata for weeks 0-52</td>
</tr>
<tr>
<td>metadata_0remission_52flare.txt</td>
<td>Metadata for remission week 0 vs flare week 52 comparison</td>
</tr>
<tr>
<td>mixOmics_drugStudy.R</td>
<td>R script for mixOmics analysis of drug study data</td>
</tr>
<tr>
<td>mixOmics_lwIBD.R</td>
<td>R script for mixOmics analysis of Living with IBD data</td>
</tr>
<tr>
<td>mixOmics.R</td>
<td>Main R script for mixOmics analysis</td>
</tr>
<tr>
<td>names</td>
<td>File containing metabolite names</td>
</tr>
<tr>
<td>parse_hmdb_xml_v2.pl</td>
<td>Perl script version 2 for parsing HMDB XML files</td>
</tr>
<tr>
<td>parse_hmdb_xml.pl</td>
<td>Perl script for parsing HMDB XML files</td>
</tr>
<tr>
<td>parse_hmdb_xml.py</td>
<td>Python script for parsing HMDB XML files</td>
</tr>
<tr>
<td>PLS-DA_SRBCT.R</td>
<td>R script for PLS-DA (Partial Least Squares Discriminant Analysis) on SRBCT dataset</td>
</tr>
<tr>
<td>PURINE_Ttest_metaboAnalyst.png</td>
<td>T-test results visualization for Purine from MetaboAnalyst</td>
</tr>
<tr>
<td>remission_samples.txt</td>
<td>List of remission sample IDs</td>
</tr>
<tr>
<td>results.tsv</td>
<td>Analysis results in TSV format</td>
</tr>
<tr>
<td>scp_dayhoff.sh</td>
<td>Shell script for SCP transfer to Dayhoff server</td>
</tr>
<tr>
<td>taxa_0remission_52flare.txt</td>
<td>Taxonomic data for remission week 0 vs flare week 52</td>
</tr>
<tr>
<td>taxa_only_0-52w_10000up_norm_remission.tsv</td>
<td>Normalized taxonomic data filtered to taxa with 10000+ counts (remission samples only)</td>
</tr>
<tr>
<td>taxa_only_0-52w_10000up_norm.tsv</td>
<td>Normalized taxonomic data filtered to taxa with 10000+ counts</td>
</tr>
<tr>
<td>taxa_only_0-52w_1000up_norm_remission.tsv</td>
<td>Normalized taxonomic data filtered to taxa with 1000+ counts (remission samples only)</td>
</tr>
<tr>
<td>taxa_only_0-52w_1000up_norm.tsv</td>
<td>Normalized taxonomic data filtered to taxa with 1000+ counts</td>
</tr>
<tr>
<td>taxa_only_0-52w_1000up.tsv</td>
<td>Taxonomic data filtered to taxa with 1000+ counts</td>
</tr>
<tr>
<td>taxa_only_0-52w_100up.tsv</td>
<td>Taxonomic data filtered to taxa with 100+ counts</td>
</tr>
<tr>
<td>taxa_only_0-52w.tsv</td>
<td>Taxonomic data for weeks 0-52</td>
</tr>
<tr>
<td>vim</td>
<td>Vim editor backup or swap file</td>
</tr>
<tr>
<td rowspan="7">0remission_52flare_analyses_230913/mixOmics/ewandson_example</td>
<td>analysis_circusplot.html</td>
<td>HTML output from circus plot analysis</td>
</tr>
<tr>
<td>analysis_circusplot.Rmd</td>
<td>R Markdown script for creating circus plots (Ewandson example)</td>
</tr>
<tr>
<td>helper_functions.R</td>
<td>Helper functions R script (Ewandson example)</td>
</tr>
<tr>
<td>results/</td>
<td>Directory containing analysis results</td>
</tr>
<tr>
<td>results1_1_1_1_cor_matrix_circos.csv</td>
<td>Correlation matrix results for circos plot</td>
</tr>
<tr>
<td>SGCCA_HFHS.csv</td>
<td>SGCCA results for HFHS dataset (Ewandson example)</td>
</tr>
<tr>
<td>SGCCA_SD.csv</td>
<td>SGCCA results for SD dataset (Ewandson example)</td>
</tr>
<tr>
<td rowspan="3">0remission_52flare_analyses_230913/mixOmics/manuscript/figures</td>
<td>20240203_LwIBD_manuscript_figures.pptx</td>
<td>PowerPoint presentation with manuscript figures from February 3, 2024</td>
</tr>
<tr>
<td rowspan="1">0remission_52flare_analyses_230913/mixOmics/manuscript/figures/juan</td>
<td>todo</td>
<td>Todo list file</td>
</tr>
<tr>
<td rowspan="3">0remission_52flare_analyses_230913/mixOmics/manuscript/figures/materials</td>
<td>correlationHeatmap_metaboAbalyst_remission0-flare52_edited.png</td>
<td>Edited correlation heatmap from MetaboAnalyst for remission week 0 vs flare week 52</td>
</tr>
<tr>
<td>correlationHeatmap_metaboAbalyst_remission0-flare52.png</td>
<td>Correlation heatmap from MetaboAnalyst for remission week 0 vs flare week 52</td>
</tr>
<tr>
<td>names.psd</td>
<td>Photoshop document with figure names/labels</td>
</tr>
<tr>
<td rowspan="4">0remission_52flare_analyses_230913/mixOmics/metaboAnalyst5.0</td>
<td>cow_diet.csv</td>
<td>Example dataset: cow diet data in CSV format</td>
</tr>
<tr>
<td>human_cachexia.csv</td>
<td>Example dataset: human cachexia data in CSV format</td>
</tr>
<tr>
<td>TCE_feature_table.csv</td>
<td>Example dataset: TCE feature table in CSV format</td>
</tr>
<tr>
<td>TCE_metadata.csv</td>
<td>Example dataset: TCE metadata in CSV format</td>
</tr>
<tr>
<td rowspan="2">0remission_52flare_analyses_230913/results_heather_231003</td>
<td>correlation_circos_remission0-flare52.png</td>
<td>Correlation circos plot for remission week 0 vs flare week 52</td>
</tr>
<tr>
<td rowspan="10">0remission_52flare_analyses_230913/results_heather_231003/livingWithIBD_Maaslin2_0remission-52flare_230927</td>
<td>all_results.tsv</td>
<td>All Maaslin2 analysis results in TSV format</td>
</tr>
<tr>
<td>diagnosis.pdf</td>
<td>PDF visualization of Maaslin2 results by diagnosis</td>
</tr>
<tr>
<td>fitted.rds</td>
<td>R data object containing fitted model from Maaslin2</td>
</tr>
<tr>
<td>maaslin2.log</td>
<td>Log file from Maaslin2 analysis</td>
</tr>
<tr>
<td>residuals.rds</td>
<td>R data object containing residuals from Maaslin2 model</td>
</tr>
<tr>
<td>significant_results.tsv</td>
<td>Significant results from Maaslin2 analysis</td>
</tr>
<tr>
<td rowspan="10">0remission_52flare_analyses_230913/results_heather_231003/livingWithIBD_Maaslin2_0remission-52flare_230927/figures</td>
<td>diagnosis_1.png</td>
<td>Diagnosis visualization figure 1</td>
</tr>
<tr>
<td>diagnosis_2.png</td>
<td>Diagnosis visualization figure 2</td>
</tr>
<tr>
<td>diagnosis_3.png</td>
<td>Diagnosis visualization figure 3</td>
</tr>
<tr>
<td>diagnosis_4.png</td>
<td>Diagnosis visualization figure 4</td>
</tr>
<tr>
<td>diagnosis_5.png</td>
<td>Diagnosis visualization figure 5</td>
</tr>
<tr>
<td>diagnosis_6.png</td>
<td>Diagnosis visualization figure 6</td>
</tr>
<tr>
<td>diagnosis_7.png</td>
<td>Diagnosis visualization figure 7</td>
</tr>
<tr>
<td>diagnosis_8.png</td>
<td>Diagnosis visualization figure 8</td>
</tr>
<tr>
<td>diagnosis_9.png</td>
<td>Diagnosis visualization figure 9</td>
</tr>
<tr>
<td>diagnosis_10.png</td>
<td>Diagnosis visualization figure 10</td>
</tr>
<tr>
<td rowspan="15">0remission_52flare_analyses_230913/results_heather_231003/livingWithIBD_Maaslin2_all0-all52_230927</td>
<td>age_group.pdf</td>
<td>PDF visualization of Maaslin2 results by age group</td>
</tr>
<tr>
<td>all_results.tsv</td>
<td>All Maaslin2 analysis results in TSV format</td>
</tr>
<tr>
<td>diagnosis.pdf</td>
<td>PDF visualization of Maaslin2 results by diagnosis</td>
</tr>
<tr>
<td>disease_severity.pdf</td>
<td>PDF visualization of Maaslin2 results by disease severity</td>
</tr>
<tr>
<td>fitted.rds</td>
<td>R data object containing fitted model from Maaslin2</td>
</tr>
<tr>
<td>heatmap.pdf</td>
<td>Heatmap visualization in PDF format</td>
</tr>
<tr>
<td>maaslin2.log</td>
<td>Log file from Maaslin2 analysis</td>
</tr>
<tr>
<td>residuals.rds</td>
<td>R data object containing residuals from Maaslin2 model</td>
</tr>
<tr>
<td>sex_f.pdf</td>
<td>PDF visualization of Maaslin2 results by sex (female)</td>
</tr>
<tr>
<td>significant_results.tsv</td>
<td>Significant results from Maaslin2 analysis</td>
</tr>
<tr>
<td>time_severity.pdf</td>
<td>PDF visualization of Maaslin2 results by time and severity</td>
</tr>
<tr>
<td>week.pdf</td>
<td>PDF visualization of Maaslin2 results by week</td>
</tr>
<tr>
<td rowspan="61">0remission_52flare_analyses_230913/results_heather_231003/livingWithIBD_Maaslin2_all0-all52_230927/figures</td>
<td>age_group_1.png</td>
<td>Age group visualization figure 1</td>
</tr>
<tr>
<td>age_group_2.png</td>
<td>Age group visualization figure 2</td>
</tr>
<tr>
<td>age_group_3.png</td>
<td>Age group visualization figure 3</td>
</tr>
<tr>
<td>age_group_4.png</td>
<td>Age group visualization figure 4</td>
</tr>
<tr>
<td>age_group_5.png</td>
<td>Age group visualization figure 5</td>
</tr>
<tr>
<td>age_group_6.png</td>
<td>Age group visualization figure 6</td>
</tr>
<tr>
<td>age_group_7.png</td>
<td>Age group visualization figure 7</td>
</tr>
<tr>
<td>age_group_8.png</td>
<td>Age group visualization figure 8</td>
</tr>
<tr>
<td>age_group_9.png</td>
<td>Age group visualization figure 9</td>
</tr>
<tr>
<td>age_group_10.png</td>
<td>Age group visualization figure 10</td>
</tr>
<tr>
<td>diagnosis_1.png</td>
<td>Diagnosis visualization figure 1</td>
</tr>
<tr>
<td>diagnosis_2.png</td>
<td>Diagnosis visualization figure 2</td>
</tr>
<tr>
<td>diagnosis_3.png</td>
<td>Diagnosis visualization figure 3</td>
</tr>
<tr>
<td>diagnosis_4.png</td>
<td>Diagnosis visualization figure 4</td>
</tr>
<tr>
<td>diagnosis_5.png</td>
<td>Diagnosis visualization figure 5</td>
</tr>
<tr>
<td>diagnosis_6.png</td>
<td>Diagnosis visualization figure 6</td>
</tr>
<tr>
<td>diagnosis_7.png</td>
<td>Diagnosis visualization figure 7</td>
</tr>
<tr>
<td>diagnosis_8.png</td>
<td>Diagnosis visualization figure 8</td>
</tr>
<tr>
<td>diagnosis_9.png</td>
<td>Diagnosis visualization figure 9</td>
</tr>
<tr>
<td>diagnosis_10.png</td>
<td>Diagnosis visualization figure 10</td>
</tr>
<tr>
<td>disease_severity_1.png</td>
<td>Disease severity visualization figure 1</td>
</tr>
<tr>
<td>disease_severity_2.png</td>
<td>Disease severity visualization figure 2</td>
</tr>
<tr>
<td>disease_severity_3.png</td>
<td>Disease severity visualization figure 3</td>
</tr>
<tr>
<td>disease_severity_4.png</td>
<td>Disease severity visualization figure 4</td>
</tr>
<tr>
<td>disease_severity_5.png</td>
<td>Disease severity visualization figure 5</td>
</tr>
<tr>
<td>disease_severity_6.png</td>
<td>Disease severity visualization figure 6</td>
</tr>
<tr>
<td>disease_severity_7.png</td>
<td>Disease severity visualization figure 7</td>
</tr>
<tr>
<td>disease_severity_8.png</td>
<td>Disease severity visualization figure 8</td>
</tr>
<tr>
<td>disease_severity_9.png</td>
<td>Disease severity visualization figure 9</td>
</tr>
<tr>
<td>disease_severity_10.png</td>
<td>Disease severity visualization figure 10</td>
</tr>
<tr>
<td>heatmap.png</td>
<td>Heatmap visualization in PNG format</td>
</tr>
<tr>
<td>sex_f_1.png</td>
<td>Sex (female) visualization figure 1</td>
</tr>
<tr>
<td>sex_f_2.png</td>
<td>Sex (female) visualization figure 2</td>
</tr>
<tr>
<td>sex_f_3.png</td>
<td>Sex (female) visualization figure 3</td>
</tr>
<tr>
<td>sex_f_4.png</td>
<td>Sex (female) visualization figure 4</td>
</tr>
<tr>
<td>sex_f_5.png</td>
<td>Sex (female) visualization figure 5</td>
</tr>
<tr>
<td>sex_f_6.png</td>
<td>Sex (female) visualization figure 6</td>
</tr>
<tr>
<td>sex_f_7.png</td>
<td>Sex (female) visualization figure 7</td>
</tr>
<tr>
<td>sex_f_8.png</td>
<td>Sex (female) visualization figure 8</td>
</tr>
<tr>
<td>sex_f_9.png</td>
<td>Sex (female) visualization figure 9</td>
</tr>
<tr>
<td>sex_f_10.png</td>
<td>Sex (female) visualization figure 10</td>
</tr>
<tr>
<td>time_severity_1.png</td>
<td>Time and severity visualization figure 1</td>
</tr>
<tr>
<td>time_severity_2.png</td>
<td>Time and severity visualization figure 2</td>
</tr>
<tr>
<td>time_severity_3.png</td>
<td>Time and severity visualization figure 3</td>
</tr>
<tr>
<td>time_severity_4.png</td>
<td>Time and severity visualization figure 4</td>
</tr>
<tr>
<td>time_severity_5.png</td>
<td>Time and severity visualization figure 5</td>
</tr>
<tr>
<td>time_severity_6.png</td>
<td>Time and severity visualization figure 6</td>
</tr>
<tr>
<td>time_severity_7.png</td>
<td>Time and severity visualization figure 7</td>
</tr>
<tr>
<td>time_severity_8.png</td>
<td>Time and severity visualization figure 8</td>
</tr>
<tr>
<td>time_severity_9.png</td>
<td>Time and severity visualization figure 9</td>
</tr>
<tr>
<td>time_severity_10.png</td>
<td>Time and severity visualization figure 10</td>
</tr>
<tr>
<td>week_1.png</td>
<td>Week visualization figure 1</td>
</tr>
<tr>
<td>week_2.png</td>
<td>Week visualization figure 2</td>
</tr>
<tr>
<td>week_3.png</td>
<td>Week visualization figure 3</td>
</tr>
<tr>
<td>week_4.png</td>
<td>Week visualization figure 4</td>
</tr>
<tr>
<td>week_5.png</td>
<td>Week visualization figure 5</td>
</tr>
<tr>
<td>week_6.png</td>
<td>Week visualization figure 6</td>
</tr>
<tr>
<td>week_7.png</td>
<td>Week visualization figure 7</td>
</tr>
<tr>
<td>week_8.png</td>
<td>Week visualization figure 8</td>
</tr>
<tr>
<td>week_9.png</td>
<td>Week visualization figure 9</td>
</tr>
<tr>
<td>week_10.png</td>
<td>Week visualization figure 10</td>
</tr>
<tr>
<td rowspan="14">0remission_52flare_analyses_230913/results_heather_231003/metaboAnalyst_remission0-flare52</td>
<td>ACETOACETATE_metaboAnalyst_remission0-flare52.png</td>
<td>Acetoacetate metabolite visualization from MetaboAnalyst</td>
</tr>
<tr>
<td>CARNOSINE_metaboAnalyst_remission0-flare52.png</td>
<td>Carnosine metabolite visualization from MetaboAnalyst</td>
</tr>
<tr>
<td>CATECHOL_metaboAnalyst_remission0-flare52.png</td>
<td>Catechol metabolite visualization from MetaboAnalyst</td>
</tr>
<tr>
<td>correlation_table_metaboAnalyst_remission0-flare52.xlsx</td>
<td>Correlation table in Excel format from MetaboAnalyst</td>
</tr>
<tr>
<td>correlationHeatmap_metaboAbalyst_remission0-flare52.png</td>
<td>Correlation heatmap from MetaboAnalyst for remission week 0 vs flare week 52</td>
</tr>
<tr>
<td>covariate_result_metaboAnalyst.xlsx</td>
<td>Covariate analysis results in Excel format from MetaboAnalyst</td>
</tr>
<tr>
<td>LInearModel_MetaboAnalyst.png</td>
<td>Linear model visualization from MetaboAnalyst</td>
</tr>
<tr>
<td>N_ACETYL_LEUCINE_metaboAnalyst_remission0-flare52.png</td>
<td>N-Acetyl Leucine metabolite visualization from MetaboAnalyst</td>
</tr>
<tr>
<td>ortho_PLSDA_metaboanalyst.png</td>
<td>Orthogonal PLS-DA plot from MetaboAnalyst</td>
</tr>
<tr>
<td>PLSDA_metaboAnalyst.png</td>
<td>PLS-DA (Partial Least Squares Discriminant Analysis) plot from MetaboAnalyst</td>
</tr>
<tr>
<td>PURINE_metaboAnalyst_remission0-flare52.png</td>
<td>Purine metabolite visualization from MetaboAnalyst</td>
</tr>
<tr>
<td>pval_correlation_table_metaboAnalyst_remission0-flare52.tsv</td>
<td>P-value correlation table in TSV format from MetaboAnalyst</td>
</tr>
<tr>
<td>top25_metabolites_correlation.png</td>
<td>Top 25 metabolites correlation visualization</td>
</tr>
<tr>
<td>VP_metaboAnalyst_remission0-flare52.png</td>
<td>Volcano plot from MetaboAnalyst for remission week 0 vs flare week 52</td>
</tr>
</tbody>
</table>

