library(data.table)
library(ggplot2)
library(grid)
library(factoextra)
library(preprocessCore)
library(gplots)
library(qvalue)
# setwd('G:/My Drive/Flu_vaccine_proteomics_project')

## compare output from different mapDIA parameter settings
expr<-as.data.frame(fread('Protein_level_original_output.txt',header=T))
expr<-as.data.frame(fread('Protein_level_mapDIA_setting_v1.txt',header=T))
expr<-as.data.frame(fread('Protein_level_mapDIA_setting_v2.txt',header=T))
expr<-as.data.frame(fread('Protein_level_mapDIA_setting_v3.txt',header=T))

expr<-expr[,c(2:(dim(expr)[2]-2))]
expr[,1:dim(expr)[2]]<-log2(expr[,1:dim(expr)[2]])

expr_pca<-prcomp(t(as.matrix(expr)),center=T,scale=T)
summary(expr_pca)

fviz_eig(expr_pca)

col_pca<-rep('Non_pool',215)
col_pca[grep('Upool',names(expr))]<-'Pool'

pdf('PCA_plot_for_protein_level_output_from_4_mapDIA_settings.pdf',width=8,height=8)
fviz_pca_ind(expr_pca,col.ind=col_pca,label='none',pointsize=2.5,alpha=.25)
dev.off()


## normalize the preprocessed expression data
protein_expr<-as.data.frame(fread('protein_level_mapDIA_setting_v1_Upool_and_rddt_spls_rmed.txt'),header=T)
protein_expr[is.na(protein_expr)]<-1
protein_expr_matrix<-as.matrix(protein_expr[,2:dim(protein_expr)[2]])
protein_expr_ColSums<-sweep(protein_expr_matrix, 2, colSums(protein_expr_matrix)/10^9, FUN="/")
protein_expr_ColSums_Quantile<-as.data.frame(round(normalize.quantiles(protein_expr_ColSums)))
protein_expr_ColSums_Quantile$Protein<-protein_expr$Protein
names(protein_expr_ColSums_Quantile)[1:(dim(protein_expr_ColSums_Quantile)[2]-1)]<-names(protein_expr)[2:dim(protein_expr)[2]]
protein_expr_ColSums_Quantile<-protein_expr_ColSums_Quantile[,c(dim(protein_expr_ColSums_Quantile)[2],1:(dim(protein_expr_ColSums_Quantile)[2]-1))]
write.table(protein_expr_ColSums_Quantile,'protein_level_mapDIA_setting_v1_Sum_Quantile.txt',quote=F,sep='\t',row.names=F)


## combine normalized data with adjusted seroconversion infomation
df1<-as.data.frame(fread('Protein_expr_Sum_Quantile_t.txt',header=T))
df2<-as.data.frame(fread('Observed_vs._predicted_seroconversion_for_UGA4_adults.csv',sep=',',header=T))
df<-merge(df1,df2[,c(2,22,24,26,28,30)])
write.table(df,'Protein_expr_Sum_Quantile_t_with_res_SC.txt',quote=F,sep='\t',row.names=F)


## log2-transform the data, row-median-center the data
expr_res<-as.data.frame(fread('Protein_expr_Sum_Quantile_t_with_res_SC.txt',header=T))
expr_res_id<-expr_res[,1]
expr_res_expr<-expr_res[,2:(dim(expr_res)[2]-5)]
expr_res_res<-expr_res[,(dim(expr_res)[2]-4):dim(expr_res)[2]]
expr_res_expr<-log2(as.matrix(expr_res_expr))
expr_res_expr_med<-apply(expr_res_expr,2,median)
expr_res_expr<-as.data.frame(sweep(expr_res_expr, 2, expr_res_expr_med, FUN="-"))
expr_res_v2<-cbind(expr_res_id,expr_res_expr,expr_res_res)
write.table(expr_res_v2,'Protein_expr_Sum_Quantile_t_log2_medCted_with_res_SC.txt',quote=F,sep='\t',row.names=F)


## correlation of top 10 PCs with the variables.  
expr<-read.table('Protein_Expr_Sum_Quantile_t_Log2_RowMedianCentered_rank_by_Seroconversion.txt',header=T,sep='\t')
protein_IDs<-expr[,1]
expr<-expr[,c(2:dim(expr)[2])]
expr_pca<-prcomp(t(as.matrix(expr)),center=T,scale=T)
summary(expr_pca)
expr_header<-read.table('Data_for_heatmap_headers_ordered_by_Seroconversion.txt',header=T)

for (x in 1:10){
  count<-0
  for (y in c(3,2,4:13)){
    corr<-as.numeric(cor.test(expr_pca$x[,x],expr_header[,y],method='s')$estimate)
    pVal<-cor.test(expr_pca$x[,x],expr_header[,y],method='s')$p.value
    print(c(names(expr_header)[y],corr,pVal))
    
    count<-count+1
    if(count%%12==0){
      print('')
    }
  }
}

## correlation of PC1 with the batches. 
PC1_batches<-read.table('PC1_of_expression_data_association_with_batch_effects.txt',header=T)
cor.test(PC1_batches$PC1,as.numeric(as.factor(PC1_batches$Column)),method='s')
cor.test(PC1_batches$PC1,as.numeric(as.factor(PC1_batches$Set)),method='s')
cor.test(PC1_batches$PC1,PC1_batches$Run_order,method='s')

lm_PC1_batches<-lm(PC1_batches$PC1~as.numeric(as.factor(PC1_batches$Column))+as.numeric(as.factor(PC1_batches$Set))+PC1_batches$Run_order)
summary(lm_PC1_batches)

## remove PC1, just for visualization (not correlated with anything, including batches). 
X = expr
mu = colMeans(X)

Xpca = prcomp(X)

Xhat = Xpca$x[,-1] %*% t(Xpca$rotation[,-1])
Xhat = scale(Xhat, center = -mu, scale = FALSE)

Xhat<-data.frame(ID=protein_IDs,as.data.frame(Xhat))
write.table(Xhat,'Protein_Expr_Sum_Quantile_t_Log2_RowMedianCentered_rank_by_Seroconversion_removing_PC1.txt',row.names = F,quote=F,sep='\t')

## heatmap of reconstructed data for all proteins in all samples. 
pdf('Heatmap_of_reconstructed_data_as_Figure1.pdf',width=8,height=8)
colors <- c(seq(-3,3,length=5000))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 4999)

heatmap.2(as.matrix(Xhat[,1:dim(Xhat)[2]]),col=my_palette,breaks=colors,trace='none',Colv=F,dendrogram = 'none', labRow = NA, main='Log2_RowMedianCentered',distfun=function(x) as.dist(1-cor(t(x),method='s')),hclustfun=function(x) hclust(x, method="complete"))
heatmap.2(as.matrix(Xhat[,1:dim(Xhat)[2]]),col=my_palette,breaks=colors,trace='none',Colv=F,dendrogram = 'none', labRow = NA, scale='row',main='Log2_RowMedianCentered_RowZscored',distfun=function(x) as.dist(1-cor(t(x),method='s')),hclustfun=function(x) hclust(x, method="complete"))
dev.off()


## expression heatmap ordered by adjusted seroconversion
expr<-read.table('Protein_Expr_Sum_Quantile_t_Log2_RowMedianCentered_rank_by_Residual_removing_PC1.txt',header=T)
pdf('Heatmap_of_protein_expression_levels_rank_by_Residual.pdf',width=8,height=8)
colors <- c(seq(-3,3,length=5000))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 4999)

heatmap.2(as.matrix(expr[,2:dim(expr)[2]]),col=my_palette,breaks=colors,trace='none',Colv=F,labRow = NA,main='Log2_RowMedianCentered',distfun=function(x) as.dist(1-cor(t(x),method='s')),hclustfun=function(x) hclust(x, method="complete"))
heatmap.2(as.matrix(expr[,2:dim(expr)[2]]),col=my_palette,breaks=colors,trace='none',Colv=F,labRow = NA,scale='column',main='Log2_RowMedianCentered_ColumnZscored',distfun=function(x) as.dist(1-cor(t(x),method='s')),hclustfun=function(x) hclust(x, method="complete"))
dev.off()

## expression heatmap ordered by raw seroconversion
expr<-read.table('Protein_Expr_Sum_Quantile_t_Log2_RowMedianCentered_rank_by_Seroconversion_removing_PC1.txt',header=T)
pdf('Heatmap_of_protein_expression_levels_rank_by_Seroconversion.pdf',width=8,height=8)
colors <- c(seq(-3,3,length=5000))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 4999)

heatmap.2(as.matrix(expr[,2:dim(expr)[2]]),col=my_palette,breaks=colors,trace='none',Colv=F,labRow = NA,main='Log2_RowMedianCentered',distfun=function(x) as.dist(1-cor(t(x),method='s')),hclustfun=function(x) hclust(x, method="complete"))
heatmap.2(as.matrix(expr[,2:dim(expr)[2]]),col=my_palette,breaks=colors,trace='none',Colv=F,labRow = NA,scale='row',main='Log2_RowMedianCentered_ColumnZscored',distfun=function(x) as.dist(1-cor(t(x),method='s')),hclustfun=function(x) hclust(x, method="complete"))
dev.off()

## demographic and other information in a heatmap for the above 2 heatmaps
heatmap_header_1<-read.table('Data_for_heatmap_headers_ordered_by_Seroconversion.txt',header=T)
heatmap_header_2<-read.table('Data_for_heatmap_headers_ordered_by_Residual.txt',header=T)
pdf('Heatmap_headers.pdf',width=8,height=4)
colors1 <- c(seq(-7.85, 10.32,length=50))
my_palette1 <- colorRampPalette(c("white", "black"))(n = 49)

colors2 <- c(seq(-3, 22,length=50))
my_palette2 <- colorRampPalette(c("white", "black"))(n = 49)

colors3 <- c(seq(0,1,length=50))
my_palette3 <- colorRampPalette(c("white", "black"))(n = 49)

colors4 <- c(seq(4.32, 6.41,length=50))
my_palette4 <- colorRampPalette(c("white", "black"))(n = 49)

colors5 <- c(seq(4.2, 5.9,length=50))
my_palette5 <- colorRampPalette(c("white", "black"))(n = 49)

colors6 <- c(seq(0, 5,length=50))
my_palette6 <- colorRampPalette(c("white", "black"))(n = 49)

colors7 <- c(seq(9, 30,length=50))
my_palette7 <- colorRampPalette(c("white", "black"))(n = 49)

heatmap.2(as.matrix(t(heatmap_header_1[,2:3])),col=my_palette1,breaks=colors1,trace='none',Rowv=F,Colv=F,labRow = NA,main='Residual')
heatmap.2(as.matrix(t(heatmap_header_1[,3:4])),col=my_palette2,breaks=colors2,trace='none',Rowv=F,Colv=F,labRow = NA,main='Seroconversion')
heatmap.2(as.matrix(t(heatmap_header_1[,c(4,7:10,12)])),col=my_palette3,breaks=colors3,trace='none',Rowv=F,Colv=F,labRow = NA,main='Categorical factors')
heatmap.2(as.matrix(t(heatmap_header_1[,5:6])),col=my_palette4,breaks=colors4,trace='none',Rowv=F,Colv=F,labRow = NA,main='Age')
heatmap.2(as.matrix(t(heatmap_header_1[,6:7])),col=my_palette5,breaks=colors5,trace='none',Rowv=F,Colv=F,labRow = NA,main='BMI')
heatmap.2(as.matrix(t(heatmap_header_1[,11:12])),col=my_palette6,breaks=colors6,trace='none',Rowv=F,Colv=F,labRow = NA,main='Month of vaccination')
heatmap.2(as.matrix(t(heatmap_header_1[,13:12])),col=my_palette7,breaks=colors7,trace='none',Rowv=F,Colv=F,labRow = NA,main='Baseline')

heatmap.2(as.matrix(t(heatmap_header_2[,2:3])),col=my_palette1,breaks=colors1,trace='none',Rowv=F,Colv=F,labRow = NA,main='Residual')
heatmap.2(as.matrix(t(heatmap_header_2[,3:4])),col=my_palette2,breaks=colors2,trace='none',Rowv=F,Colv=F,labRow = NA,main='Seroconversion')
heatmap.2(as.matrix(t(heatmap_header_2[,c(4,7:10,12)])),col=my_palette3,breaks=colors3,trace='none',Rowv=F,Colv=F,labRow = NA,main='Categorical factors')
heatmap.2(as.matrix(t(heatmap_header_2[,5:6])),col=my_palette4,breaks=colors4,trace='none',Rowv=F,Colv=F,labRow = NA,main='Age')
heatmap.2(as.matrix(t(heatmap_header_2[,6:7])),col=my_palette5,breaks=colors5,trace='none',Rowv=F,Colv=F,labRow = NA,main='BMI')
heatmap.2(as.matrix(t(heatmap_header_2[,11:12])),col=my_palette6,breaks=colors6,trace='none',Rowv=F,Colv=F,labRow = NA,main='Month of vaccination')
heatmap.2(as.matrix(t(heatmap_header_2[,13:12])),col=my_palette7,breaks=colors7,trace='none',Rowv=F,Colv=F,labRow = NA,main='Baseline')
dev.off()


## vocano plots for differential expression based on raw & adjusted serocoversion
seroconversion_comparison<-read.table('FoldChange_pVal_Comparing_Seroconversion.txt',header=T,sep='\t')
residual_comparison<-read.table('FoldChange_pVal_Comparing_Residual.txt',header=T,sep='\t')
seroconversion_comparison_cholesterol<-seroconversion_comparison[c(2,1,96,9,8,3,17,7,5),]
residual_comparison_actin<-residual_comparison[c(56,1,3,18),]
seroconversion_comparison_cholesterol$GeneName<-c('APOA1','APOC2','LPA','APOE','APOA4','APOA2','CETP','APOC3','LCAT')
residual_comparison_actin$GeneName<-c('TMSB4X','CFL1','PFN1','ACTB')

pdf('FoldChange_pVal_for_Residual_and_Seroconversion_comparison.pdf',width=6,height=6)
p1<-ggplot(seroconversion_comparison,aes(x=log2.Fold.change.,y=-log10(pVal)))+geom_point(size=3,alpha=.5)+
  # geom_hline(yintercept = -log10(0.01),linetype=5)+
  labs(x='Fold change high/low SC (log base 2)',y='P-value (-log base 10)')+
  scale_x_continuous(limits=c(-1,1))+
  scale_y_continuous(limits=c(0,4))+
  theme_bw()+theme(text=element_text(size=15),panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

p2<-ggplot(residual_comparison,aes(x=log2.Fold.change.,y=-log10(pVal)))+geom_point(size=3,alpha=.5)+
  # geom_hline(yintercept = -log10(0.01),linetype=5)+
  labs(x='Fold change high/low SC* (log base 2)',y='P-value (-log base 10)')+
  scale_x_continuous(limits=c(-1,1))+
  scale_y_continuous(limits=c(0,4))+
  theme_bw()+theme(text=element_text(size=15),panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

p3<-ggplot(seroconversion_comparison,aes(x=log2.Fold.change.,y=-log10(pVal)))+geom_point(size=3,alpha=.5)+
  geom_point(data=seroconversion_comparison_cholesterol,aes(x=log2.Fold.change.,y=-log10(pVal)),col='red',size=5,alpha=.5)+geom_text(data=seroconversion_comparison_cholesterol,aes(label=GeneName),hjust=0.5,vjust=-0.75,col='red')+
  # geom_hline(yintercept = -log10(0.01),linetype=5)+
  labs(x='Fold change high/low SC (log base 2)',y='P-value (-log base 10)')+
  scale_x_continuous(limits=c(-1,1))+
  scale_y_continuous(limits=c(0,4))+
  theme_bw()+theme(text=element_text(size=15),panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

p4<-ggplot(residual_comparison,aes(x=log2.Fold.change.,y=-log10(pVal)))+geom_point(size=3,alpha=.5)+
  geom_point(data=residual_comparison_actin,aes(x=log2.Fold.change.,y=-log10(pVal)),col='red',size=5,alpha=.5)+geom_text(data=residual_comparison_actin,aes(label=GeneName),hjust=0.5,vjust=-0.75,col='red')+
  # geom_hline(yintercept = -log10(0.01),linetype=5)+
  labs(x='Fold change high/low SC* (log base 2)',y='P-value (-log base 10)')+
  scale_x_continuous(limits=c(-1,1))+
  scale_y_continuous(limits=c(0,4))+
  theme_bw()+theme(text=element_text(size=15),panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

print(p1)
print(p2)
print(p3)
print(p4)
dev.off()


## expression heatmaps for leading genes in 2 enriched pathways
expr_seroconversion<-read.table('Protein_Expr_Sum_Quantile_t_Log2_RowMedianCentered_rank_by_Seroconversion_removing_PC1.txt',header=T)
expr_residual<-read.table('Protein_Expr_Sum_Quantile_t_Log2_RowMedianCentered_rank_by_Residual_removing_PC1.txt',header=T)
seroconversion_sig<-read.table('Leading_genes_in_Cholesterol_Metablism_pathway.txt',header=T)
residual_sig<-read.table('Leading_genes_in_Regulation_of_Cytoskeleton_pathway.txt',header=T)
expr_seroconversion_sig<-merge(seroconversion_sig,expr_seroconversion,by='ID')
expr_residual_sig<-merge(residual_sig,expr_residual,by='ID')
rownames(expr_seroconversion_sig)<-expr_seroconversion_sig[,2]
rownames(expr_residual_sig)<-expr_residual_sig[,2]
expr_seroconversion_sig<-expr_seroconversion_sig[,-2]
expr_residual_sig<-expr_residual_sig[,-2]

pdf('Heatmap_for_2_significant_pathways.pdf',width=8,height=8)
colors <- c(seq(-3,3,length=1000))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 999)

heatmap.2(as.matrix(expr_seroconversion_sig[,c(2:31,(dim(expr_seroconversion_sig)[2]-40):(dim(expr_residual_sig)[2]-11))]),col=my_palette,breaks=colors,trace='none',Colv=F,dendrogram='row',scale='row',main='Leading genes in Cholesterol Metabolism pathway',distfun=function(x) as.dist(1-cor(t(x),method='s')),hclustfun=function(x) hclust(x, method="complete"),labCol=NA)
heatmap.2(as.matrix(expr_residual_sig[,c(2:31,(dim(expr_residual_sig)[2]-29):dim(expr_residual_sig)[2])]),col=my_palette,breaks=colors,trace='none',Colv=F,dendrogram='row',scale='row',main='Leading genes in Regulation of Cytoskeleton pathway',distfun=function(x) as.dist(1-cor(t(x),method='s')),hclustfun=function(x) hclust(x, method="complete"),labCol=NA)
dev.off()


## two-way ANOVA involving age/BMI and adjusted seroconversion.
flu_proteomics<-read.table('Age_and_BMI_specific_markers_ANOVA_dataset.txt',header=T)
for (i in 4:dim(flu_proteomics)[2]){
  inter_age<-as.numeric(unlist(summary(aov(flu_proteomics[1:50,i]~flu_proteomics[1:50,1]*flu_proteomics[1:50,2])))[19])
  inter_BMI<-as.numeric(unlist(summary(aov(flu_proteomics[51:110,i]~flu_proteomics[51:110,1]*flu_proteomics[51:110,3])))[19])
  res<-t(c(names(flu_proteomics)[i],inter_age,inter_BMI))
  write.table(res,'Age_and_BMI_specific_markers_ANOVA_results.txt',quote=F,sep='\t',row.names = F,col.names = F,append=T)
}

pval<-read.table('Age_and_BMI_specific_markers_ANOVA_results.txt',header=T)
pval[,4]<-p.adjust(pval[,2],method = 'BH')
pval[,5]<-p.adjust(pval[,3],method = 'BH')
names(pval)[4:5]<-c('qVal_Interaction_response_and_age','qVal_Interaction_response_and_BMI')
write.table(pval,'Age_and_BMI_specific_markers_ANOVA_results.txt',quote=F,sep='\t',row.names = F)

## GSEA plot and heatmaps for leading genes in the enriched pathways
gsea_results<-read.table('GSEA_enrichment_results_for_age_specific_markers.txt',header=T,sep='\t')
age_specific_expr<-read.table('Expression_data_frame_for_age_specific_markers_removing_PC1.txt',header=T)
gsea_results<-gsea_results[order(gsea_results$normalizedEnrichmentScore,decreasing = T),]
gsea_results$description<-factor(gsea_results$description,levels=rev(gsea_results$description))
gsea_results_IDs<-unique(unlist(strsplit(gsea_results$userId,';')))
gsea_results_expr<-merge(data.frame(ID=gsea_results_IDs,Rown=1:length(gsea_results_IDs)),age_specific_expr)
gsea_results_expr<-gsea_results_expr[order(gsea_results_expr$Rown,decreasing = F),]
gsea_results_expr<-gsea_results_expr[,-2]

pdf('GSEA_results_for_age_specific_markers.pdf',width=8,height=8)
p1<-ggplot(gsea_results,aes(x=normalizedEnrichmentScore,y=description))+geom_point(aes(color=FDR,size=size))+
  labs(x='Enrichment score',y=NULL)+
  theme_classic()+theme(legend.position = 'top')
print(p1)

colors <- c(seq(-3,3,length=300))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

heatmap.2(as.matrix(gsea_results_expr[1:26,2:dim(gsea_results_expr)[2]]),col=my_palette,breaks=colors,trace='none',Colv=F,dendrogram='row',scale='row',main='Leading genes in GSEA enriched pathways',distfun=function(x) as.dist(1-cor(t(x),method='s')),hclustfun=function(x) hclust(x, method="complete"),labCol=NA,labRow=NA)
heatmap.2(as.matrix(gsea_results_expr[27:60,2:dim(gsea_results_expr)[2]]),col=my_palette,breaks=colors,trace='none',Colv=F,dendrogram='row',scale='row',main='Leading genes in GSEA enriched pathways',distfun=function(x) as.dist(1-cor(t(x),method='s')),hclustfun=function(x) hclust(x, method="complete"),labCol=NA,labRow=NA)
dev.off()


## strain-specific signatures of the vaccine response
library(UpSetR)
sig_proteins <- read.table('Significant_proteins_in_individual_strains.txt', header=T)
all_proteins <- read.table('FoldChange_pVal_Comparing_Residual_for_individual_strains.txt',header=T,sep='\t')
all_proteins_H1N1_unique<-all_proteins[all_proteins$pVal_H1N1<=.05 & all_proteins$pVal_H3N2>.05 & all_proteins$pVal_IBV_Yam>.05 & all_proteins$pVal_IBV_Vic>.05,]
all_proteins_H3N2_unique<-all_proteins[all_proteins$pVal_H1N1>.05 & all_proteins$pVal_H3N2<=.05 & all_proteins$pVal_IBV_Yam>.05 & all_proteins$pVal_IBV_Vic>.05,]
all_proteins_IBV_Yam_unique<-all_proteins[all_proteins$pVal_H1N1>.05 & all_proteins$pVal_H3N2>.05 & all_proteins$pVal_IBV_Yam<=.05 & all_proteins$pVal_IBV_Vic>.05,]
all_proteins_IBV_Vic_unique<-all_proteins[all_proteins$pVal_H1N1>.05 & all_proteins$pVal_H3N2>.05 & all_proteins$pVal_IBV_Yam>.05 & all_proteins$pVal_IBV_Vic<=.05,]
all_proteins_common_H1N1_IBV_Vic<-all_proteins[all_proteins$pVal_H1N1<=.05 & all_proteins$pVal_H3N2>.05 & all_proteins$pVal_IBV_Yam>.05 & all_proteins$pVal_IBV_Vic<=.05,]

pdf('Strain_specific_proteomic_expressions.pdf',width=6,height=6)
p1<-upset(sig_proteins[,c(1,rev(2:11))], sets = rev(names(sig_proteins)[2:11]), sets.bar.color = "#56B4E9", keep.order = T)

p2<-ggplot(all_proteins,aes(x=log2..Fold.change._H1N1,y=-log10(pVal_H1N1)))+geom_point(size=3,alpha=.5)+
  geom_point(data=all_proteins[all_proteins$pVal_H1N1<=.05,],aes(x=log2..Fold.change._H1N1,y=-log10(pVal_H1N1)),size=5,alpha=.5)+geom_text(data=all_proteins[all_proteins$pVal_H1N1<=.05,],aes(label=GeneSymbol),hjust=0.5,vjust=-0.75)+
  geom_point(data=all_proteins_H1N1_unique,aes(x=log2..Fold.change._H1N1,y=-log10(pVal_H1N1)),col='red',size=5,alpha=.5)+geom_text(data=all_proteins_H1N1_unique,aes(label=GeneSymbol),hjust=0.5,vjust=-0.75,col='red')+
  geom_point(data=all_proteins_common_H1N1_IBV_Vic,aes(x=log2..Fold.change._H1N1,y=-log10(pVal_H1N1)),col='blue',size=5,alpha=.5)+geom_text(data=all_proteins_common_H1N1_IBV_Vic,aes(label=GeneSymbol),hjust=0.5,vjust=-0.75,col='blue')+
  # geom_hline(yintercept = -log10(0.01),linetype=5)+
  labs(x='Fold change high/low SC* (log base 2)',y='P-value (-log base 10)',title='H1N1')+
  scale_x_continuous(limits=c(-1,1.05))+
  scale_y_continuous(limits=c(0,2.75))+
  theme_bw()+theme(text=element_text(size=15),panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

p3<-ggplot(all_proteins,aes(x=log2..Fold.change._H3N2,y=-log10(pVal_H3N2)))+geom_point(size=3,alpha=.5)+
  geom_point(data=all_proteins[all_proteins$pVal_H3N2<=.05,],aes(x=log2..Fold.change._H3N2,y=-log10(pVal_H3N2)),size=5,alpha=.5)+geom_text(data=all_proteins[all_proteins$pVal_H3N2<=.05,],aes(label=GeneSymbol),hjust=0.5,vjust=-0.75)+
  geom_point(data=all_proteins_H3N2_unique,aes(x=log2..Fold.change._H3N2,y=-log10(pVal_H3N2)),col='red',size=5,alpha=.5)+geom_text(data=all_proteins_H3N2_unique,aes(label=GeneSymbol),hjust=0.5,vjust=-0.75,col='red')+
  # geom_hline(yintercept = -log10(0.01),linetype=5)+
  labs(x='Fold change high/low SC* (log base 2)',y='P-value (-log base 10)',title='H3N2')+
  scale_x_continuous(limits=c(-1,1.05))+
  scale_y_continuous(limits=c(0,2.75))+
  theme_bw()+theme(text=element_text(size=15),panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

p4<-ggplot(all_proteins,aes(x=log2..Fold.change._IBV_Yam,y=-log10(pVal_IBV_Yam)))+geom_point(size=3,alpha=.5)+
  geom_point(data=all_proteins[all_proteins$pVal_IBV_Yam<=.05,],aes(x=log2..Fold.change._IBV_Yam,y=-log10(pVal_IBV_Yam)),size=5,alpha=.5)+geom_text(data=all_proteins[all_proteins$pVal_IBV_Yam<=.05,],aes(label=GeneSymbol),hjust=0.5,vjust=-0.75)+
  geom_point(data=all_proteins_IBV_Yam_unique,aes(x=log2..Fold.change._IBV_Yam,y=-log10(pVal_IBV_Yam)),col='red',size=5,alpha=.5)+geom_text(data=all_proteins_IBV_Yam_unique,aes(label=GeneSymbol),hjust=0.5,vjust=-0.75,col='red')+
  # geom_hline(yintercept = -log10(0.01),linetype=5)+
  labs(x='Fold change high/low SC* (log base 2)',y='P-value (-log base 10)',title='IBV_Yam')+
  scale_x_continuous(limits=c(-1,1.05))+
  scale_y_continuous(limits=c(0,2.75))+
  theme_bw()+theme(text=element_text(size=15),panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

p5<-ggplot(all_proteins,aes(x=log2..Fold.change._IBV_Vic,y=-log10(pVal_IBV_Vic)))+geom_point(size=3,alpha=.5)+
  geom_point(data=all_proteins[all_proteins$pVal_IBV_Vic<=.05,],aes(x=log2..Fold.change._IBV_Vic,y=-log10(pVal_IBV_Vic)),size=5,alpha=.5)+geom_text(data=all_proteins[all_proteins$pVal_IBV_Vic<=.05,],aes(label=GeneSymbol),hjust=0.5,vjust=-0.75)+
  geom_point(data=all_proteins_IBV_Vic_unique,aes(x=log2..Fold.change._IBV_Vic,y=-log10(pVal_IBV_Vic)),col='red',size=5,alpha=.5)+geom_text(data=all_proteins_IBV_Vic_unique,aes(label=GeneSymbol),hjust=0.5,vjust=-0.75,col='red')+
  geom_point(data=all_proteins_common_H1N1_IBV_Vic,aes(x=log2..Fold.change._IBV_Vic,y=-log10(pVal_IBV_Vic)),col='blue',size=5,alpha=.5)+geom_text(data=all_proteins_common_H1N1_IBV_Vic,aes(label=GeneSymbol),hjust=0.5,vjust=-0.75,col='blue')+
  # geom_hline(yintercept = -log10(0.01),linetype=5)+
  labs(x='Fold change high/low SC* (log base 2)',y='P-value (-log base 10)',title='IBV_Vic')+
  scale_x_continuous(limits=c(-1,1.05))+
  scale_y_continuous(limits=c(0,2.75))+
  theme_bw()+theme(text=element_text(size=15),panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()


## protein expression in lean and obese for genes LRG and FHR-4. 
expr<-read.table('LRG_and_FHR4_expression_in_lean_and_obese.txt',header=T)
expr_t<-as.data.frame(t(expr[,2:dim(expr)[2]]))
names(expr_t)<-expr[,1]
expr_t$Group<-c(rep('Lean',31),rep('Obese',63))
t.test(expr_t$P02750[expr_t$Group=='Lean'],expr_t$P02750[expr_t$Group=='Obese']) ## LRG, P = 0.004. 
t.test(expr_t$Q92496[expr_t$Group=='Lean'],expr_t$Q92496[expr_t$Group=='Obese']) ## FHR4, P = 0.03. 

pdf('Expression_of_LRG_and_FHR4_in_lean_and_obese_groups.pdf',width=6,height=8)
p1<-ggplot(expr_t,aes(x=Group,y=P02750))+geom_boxplot(width=.6)+
  labs(y='LRG expression')+
  theme_classic()+theme(axis.title.x = element_blank())
p2<-ggplot(expr_t,aes(x=Group,y=Q92496))+geom_boxplot(width=.6)+
  labs(y='FHR4 expression')+
  theme_classic()+theme(axis.title.x = element_blank())
print(p1)
print(p2)
dev.off()


## correlation-based analysis. 
sc<-read.table('Data_for_heatmap_headers_ordered_by_Seroconversion.txt',header=T)
expr<-read.table('Protein_Expr_Sum_Quantile_t_Log2_RowMedianCentered_rank_by_Seroconversion.txt',header=T)
for (x in 1:dim(expr)[1]){
  sc1_expr_corr<-cor.test(sc[,3],t(expr[x,2:dim(expr)[2]]),method='s')
  sc2_expr_corr<-cor.test(sc[,2],t(expr[x,2:dim(expr)[2]]),method='s')
  
  sc2_expr_lm<-lm(sc[,2]~t(expr[x,2:dim(expr)[2]])+sc[,5]+t(expr[x,2:dim(expr)[2]])*sc[,5])
  
  age_expr_corr<-cor.test(sc[,5],t(expr[x,2:dim(expr)[2]]),method='s')
  bmi_expr_corr<-cor.test(sc[,6],t(expr[x,2:dim(expr)[2]]),method='s')
  expr_age_bmi_lm<-lm(t(expr[x,2:dim(expr)[2]])~sc[,5]+sc[,6]+sc[,7]+sc[,5]*sc[,6])
  
  result<-t(c(expr$ID[x],sc1_expr_corr$estimate,sc1_expr_corr$p.value,sc2_expr_corr$estimate,sc2_expr_corr$p.value,as.vector(t(summary(sc2_expr_lm)$coefficients[2:4,c(1,4)])),age_expr_corr$estimate,age_expr_corr$p.value,bmi_expr_corr$estimate,bmi_expr_corr$p.value,as.vector(t(summary(expr_age_bmi_lm)$coefficients[2:5,c(1,4)]))))
  write.table(result,'Correlation_based_analysis_results.txt',quote=F,sep='\t',row.names=F,col.names=F,append=T)
}

## plot the results. 
corr<-read.table('Correlation_based_analysis_results.txt',header=T)
cholesterol_pathway_ids<-read.table('Leading_genes_in_Cholesterol_Metablism_pathway.txt',header=T)
actin_pathway_ids<-read.table('Leading_genes_in_Regulation_of_Cytoskeleton_pathway.txt',header=T)
corr_cholesterol<-merge(corr,cholesterol_pathway_ids,by='ID')
corr_actin<-merge(corr,actin_pathway_ids,by='ID')

pdf('Correlation_analysis_results_plots.pdf',width=6,height=6)
p1<-ggplot(corr,aes(x=SC_Expr_Rs_coef,y=-log10(SC_Expr_Rs_pval)))+geom_point(size=3,alpha=.5)+
  geom_point(data=corr_cholesterol,aes(x=SC_Expr_Rs_coef,y=-log10(SC_Expr_Rs_pval)),col='red',size=5,alpha=.5)+geom_text(data=corr_cholesterol,aes(label=Symbol),hjust=0.5,vjust=-0.75,col='red')+
  labs(x='Corr. coef. (Seroconversion vs. Expr.)',y='P-value (-log base 10)')+
  scale_x_continuous(limits=c(-0.35,0.35))+
  scale_y_continuous(limits=c(0,4.5))+
  theme_bw()+theme(text=element_text(size=15),panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

p2<-ggplot(corr,aes(x=SC._Expr_Rs_coef,y=-log10(SC._Expr_Rs_pval)))+geom_point(size=3,alpha=.5)+
  geom_point(data=corr_actin,aes(x=SC._Expr_Rs_coef,y=-log10(SC._Expr_Rs_pval)),col='red',size=5,alpha=.5)+geom_text(data=corr_actin,aes(label=Symbol),hjust=0.5,vjust=-0.75,col='red')+
  labs(x='Corr. coef. (Seroconversion* vs. Expr.)',y='P-value (-log base 10)')+
  scale_x_continuous(limits=c(-0.35,0.35))+
  scale_y_continuous(limits=c(0,4.5))+
  theme_bw()+theme(text=element_text(size=15),panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())
print(p1)
print(p2)
dev.off()






