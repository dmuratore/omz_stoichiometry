## Code for generating stoichiometric profiles and 
## calculating ordinations and post hoc models
## for main text Figures 2 and 3
## Gene-level stoichiometric characteristics were calculated
## using the python script stoich_characteristics.py, based off
## of code provided by Mende et al 2017
## https://github.com/JessAwBryant/gene-characteristics
## Gene annotations conducted using GhostKOALA KEGG annotation server
## Coverages and estimated copy numbers calculated
## using working_with_coverage.R 
## Copy numbers estimated by mean ratio to 
## coverage the set of COGs "COG0012","COG0016",
## "COG0018","COG0172","COG0215","COG0495",
## "COG0525","COG0533","COG0541","COG0552"


## Loading libraries
library(data.table)
library(vegan)
library(tidyverse)
library(patchwork)
library(nlme)

## Read in gene-by-gene coverage information
all_gene_coverage<-fread('data/gene_mapping/full_cov_table.csv',header=TRUE,sep=',')

## Now we get gene stoichiometric data w/eggNOG functional annotations
gene_data<-fread('data/gene_mapping/full_gene_info.csv',drop=c(1,2))

## Separating out the estimated per-genome copy numbers from the gene coverages
just_cn<-all_gene_coverage[,c(seq(2,116,by=2))]
just_cov<-all_gene_coverage[,c(seq(1,115,by=2))]

## Getting row names back on to gene-by-gene coverages
rownames(all_gene_coverage)<-gene_data$gene.id
rownames(just_cn)<-gene_data$gene.id
rownames(just_cov)<-gene_data$gene.id

## Getting the total per-metagenome coverages in reads/bp gene length
coverage_sums<-rowSums(just_cov)

## Now setting up how to get the subsets we want
get_subsets<-function(df){
  subset_coverages<-just_cov[which(rownames(just_cov) %in% df$gene.id),]
  sum_coverages<-colSums(subset_coverages) ## Getting total coverages summed over all genes
  scaled_coverages<-scale(subset_coverages) ## Scaling for intercomparability
  weighted_cov_mat<-t(t(subset_coverages)/colSums(subset_coverages)) ## Relative coverage of each viral gene per bp covered
  stat_mat<-as.matrix(df[,c(2:29)])
  results<-list()
  for(i in 1:ncol(stat_mat)){
    temp<-as.numeric(stat_mat[,i])*weighted_cov_mat
    results[[i]]<-colSums(temp,na.rm=TRUE)
  }
  summary_mat<-do.call(cbind,results)
  colnames(summary_mat)<-colnames(stat_mat)
  samp_names<-rownames(summary_mat)
  year<-gsub('ETNP','',samp_names)
  year<-gsub('S.*','',year)
  fractions<-rep(NA,length(year))
  fractions[grep('SV',samp_names)]<-'SV'
  fractions[grep('PF',samp_names)]<-'PF'
  deeps<-gsub('^.*MG','',samp_names)
  deeps<-gsub('[A-Z].*$','',deeps)
  final_summary_mat<-data.frame(summary_mat,
                                year,fractions,deeps)
  return(list(df,subset_coverages,sum_coverages,scaled_coverages,weighted_cov_mat,final_summary_mat))
}

## Separate genes annotated as viral, archaeal, bacterial using the above function

virus_information<-get_subsets(filter(gene_data,broad_tax=='Viruses'))
arc_information<-get_subsets(filter(gene_data,broad_tax=='Archaea'))
bac_information<-get_subsets(filter(gene_data,broad_tax=='Bacteria'))

## Function for conducting ordination, 
## requires reading files from omz_stoich_f1_script.R
make_ordination_ready<-function(list){
  final_summary_mat<-list[[6]]
  # Parsing the depths to be modeled
  final_summary_mat$deeps<-as.numeric(as.character(final_summary_mat$deeps))
  final_summary_mat <- final_summary_mat %>%
    mutate(depth_class=round(deeps,-1))
  # Preparing stoichiogenomic parameter matrix for ordination
  ordination_mat<-t(final_summary_mat)
  ordination_mat<-ordination_mat[1:28,]
  ordination_mat<-apply(ordination_mat,1,as.numeric)
  # Scaling
  ordination_mat_scale<-t(scale(ordination_mat))
  colnames(ordination_mat_scale)<-rownames(final_summary_mat)
  
  ## Purpose here is to remove the confounding of sequencing depth with other variables
  cca_ord<-vegan::rda(t(ordination_mat_scale)~nreads,data=bulk_mat[which(bulk_mat$sample %in% gsub('_cov','',colnames(ordination_mat_scale))),])
  ## Post hoc tests using size fraction and depth to investigate relative explanatory power
  adonis2(t(ordination_mat_scale)~nreads+fraction+asinh(depth),
          data=bulk_mat[which(bulk_mat$sample %in% gsub('_cov',
                                                        '',
                                                        colnames(ordination_mat_scale))),],
          method='euclidean')
  ## Preparing outputs to be saved for figure
  ord_cols<-vegan::scores(cca_ord,choices=c(2,3))$sites
  colnames(ord_cols)<-c(summary(cca_ord)$cont$importance[2,2:3])
  ca_frame<-cbind(final_summary_mat,ord_cols)
  return(ca_frame)
}

## Conducting analyses
## Bacterial genes
bac_pca<-make_ordination_ready(bac_information)
pc_mags<-round(as.numeric(colnames(bac_pca)[33:34]),4)*100
colnames(bac_pca)[33:34]<-c('PC1','PC2')

b_fig<-ggplot(bac_pca,aes(x=PC1,y=PC2,fill=deeps,shape=fractions))+
  geom_point(size=3,col='black')+
  theme_bw()+
  scale_fill_gradient2(low='chocolate2',mid='white',high='darkblue',
                       midpoint=log10(80),trans='log10',name='Depth [m]')+
  scale_shape_manual(values=c(21,24),
                     name='Size Fraction',
                     labels=c('Particle','Free-Living'),
                     guide=guide_legend(nrow=2))+
  coord_fixed(ratio=pc_mags[2]/pc_mags[1])+
  ylab(paste0('PC1 ',pc_mags[2],'%'))+
  xlab(paste0('PC2 ',pc_mags[1],'%'))+
  ggtitle('Bacterial Genes')+
  theme(legend.position='none')

## Archaeal genes
arc_pca<-make_ordination_ready(arc_information)
pc_mags<-round(as.numeric(colnames(arc_pca)[33:34]),4)*100
colnames(arc_pca)[33:34]<-c('PC1','PC2')

a_fig<-ggplot(arc_pca,aes(x=PC1,y=PC2,fill=deeps,shape=fractions))+
  geom_point(size=3,col='black')+
  theme_bw()+
  scale_fill_gradient2(low='chocolate2',mid='white',high='darkblue',
                       midpoint=log10(80),trans='log10',name='Depth [m]')+
  scale_shape_manual(values=c(21,24),
                     name='Size Fraction',
                     labels=c('Particle','Free-Living'),
                     guide=guide_legend(nrow=2))+
  coord_fixed(ratio=pc_mags[2]/pc_mags[1])+
  ylab(paste0('PC1 ',pc_mags[2],'%'))+
  xlab(paste0('PC2 ',pc_mags[1],'%'))+
  ggtitle('Archaeal Genes')+
  theme(legend.position='none')

## Viral genes
viral_pca<-make_ordination_ready(virus_information)
pc_mags<-round(as.numeric(colnames(viral_pca)[33:34]),4)*100
colnames(viral_pca)[33:34]<-c('PC1','PC2')

v_fig<-ggplot(viral_pca,aes(x=PC1,y=PC2,fill=deeps,shape=fractions))+
  geom_point(size=3,col='black')+
  theme_bw()+
  scale_fill_gradient2(low='chocolate2',mid='white',high='darkblue',
                       midpoint=log10(80),trans='log10',name='Depth [m]')+
  scale_shape_manual(values=c(21,24),
                     name='Size Fraction',
                     labels=c('Particle','Free-Living'),
                     guide=guide_legend(nrow=2))+
  coord_fixed(ratio=pc_mags[2]/pc_mags[1])+
  ylab(paste0('PC1 ',pc_mags[2],'%'))+
  xlab(paste0('PC2 ',pc_mags[1],'%'))+
  theme(legend.position='bottom')+
  ggtitle('Viral Genes')

## Composing final figure
pca_fig<-rbind(ggplotGrob(b_fig),ggplotGrob(a_fig),ggplotGrob(v_fig),size='last')
ggsave('figures/F2.pdf',pca_fig,device='pdf')


gene_type<-c(rep('Bacteria',nrow(bac_pca)),
             rep('Archaea',nrow(arc_pca)),
             rep('Virus',nrow(viral_pca)))
all_domains_frame<-data.frame(rbind(bac_pca,arc_pca,viral_pca),gene_type)

## Estimating genome lengths based on relative abundance of single
## copy core genes
## Using COG0012 as case study gene
cog_fasta<-fread('data/sequences/cog12.fna',sep='>',header=FALSE)
cog_ids<-cog_fasta[grep('>',cog_fasta$V1),]
cog_id_nocarat<-gsub(' ','',gsub('>','',cog_ids$V1))

## Get coverages from metagenomes
cog_info<-bac_information[[1]][which(bac_information[[1]]$gene.id %in% cog_id_nocarat),]

## Get average copy number from each sample
w_cpn<-just_cn %>%
  rownames_to_column() %>%
  right_join(.,bac_information[[1]],by=c('rowname'='gene.id'))

## Matrix multiplication to get average gene length scaled by estimated
## COG copy number summed together
multiplied_out<-t(w_cpn$gene_length) %*% as.matrix(as.data.frame(w_cpn)[,grep('cn',colnames(w_cpn))]) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(adj_rowname=gsub('_cn','',rowname))

## Adding on to bacterial PCA data
new_pca_frame<-data.frame(bac_pca,est_g_length=multiplied_out$V1)

## Generating estimated genome length box plot
egl_box<-ggplot(new_pca_frame)+
  theme_bw()+
  geom_boxplot(aes(x=fractions,y=est_g_length/1e6))+
  geom_jitter(aes(x=fractions,y=est_g_length/1e6,fill=deeps),
              shape=21,
              col='black',
              width=0.25,
              size=3)+
  scale_fill_gradient2(low='chocolate2',mid='white',high='darkblue',
                       midpoint=log10(80),trans='log10',name='Depth [m]')+
  scale_x_discrete(labels=c('Particle','Free\nliving'))+
  guides(fill=guide_colorbar(reverse=TRUE))+
  ylab('Estimated Genome Length Mbp')+
  theme(strip.background=element_blank(),
        text=element_text(size=16),
        axis.title.x=element_blank())

## Plotting stoichiogenomic parameters for Figure 3
gc_box<-ggplot(all_domains_frame,aes(x=fractions))+
  theme_bw()+
  geom_boxplot(aes(y=gene_GC))+
  geom_jitter(aes(y=gene_GC,fill=deeps),
              shape=21,
              col='black',
              width=0.25,
              size=3)+
  facet_wrap(~gene_type)+
  scale_fill_gradient2(low='chocolate2',mid='white',high='darkblue',
                       midpoint=log10(80),trans='log10',name='Depth [m]')+
  scale_x_discrete(labels=c('Particle','Free\nliving'))+
  guides(fill=guide_colorbar(reverse=TRUE))+
  ylab('GC %')+
  theme(strip.background=element_blank(),
        text=element_text(size=16),
        axis.title.x=element_blank(),
        legend.position='none')
narsc_box<-ggplot(all_domains_frame,aes(x=fractions))+
  theme_bw()+
  geom_boxplot(aes(y=N_ARSC))+
  geom_jitter(aes(y=N_ARSC,fill=deeps),
              shape=21,
              col='black',
              width=0.25,
              size=3)+
  facet_wrap(~gene_type)+
  scale_fill_gradient2(low='chocolate2',mid='white',high='darkblue',
                       midpoint=log10(80),trans='log10',name='Depth [m]')+
  scale_x_discrete(labels=c('Particle','Free\nliving'))+
  guides(fill=guide_colorbar(reverse=TRUE))+
  ylab('Avg #N per AA Side Chain')+
  theme(strip.background=element_blank(),
        text=element_text(size=16),
        axis.title.x=element_blank(),
        legend.position='none')
nc_box<-ggplot(all_domains_frame,aes(x=fractions))+
  theme_bw()+
  geom_boxplot(aes(y=N_C_rat))+
  geom_jitter(aes(y=N_C_rat,fill=deeps),
              shape=21,
              col='black',
              width=0.25,
              size=3)+
  facet_wrap(~gene_type)+
  scale_fill_gradient2(low='chocolate2',mid='white',high='darkblue',
                       midpoint=log10(80),trans='log10',name='Depth [m]')+
  scale_x_discrete(labels=c('Particle','Free\nliving'))+
  guides(fill=guide_colorbar(reverse=TRUE))+
  ylab('Protein N:C Ratio')+
  theme(strip.background=element_blank(),
        text=element_text(size=16),
        axis.title.x=element_blank(),
        legend.position='none')

all_box<-(gc_box/narsc_box/nc_box)|egl_box

ggsave('figures/F3.pdf',device='pdf',
       all_box,scale=1.5)

## Conducting mixed model analysis

## For bacterial genome length between size fractions
for_modeling<-data.frame(new_pca_frame,nreads=bulk_mat$nreads)
length_model<-lm(data=for_modeling,formula=est_g_length~nreads+fractions)
wilcox.test(est_g_length~fractions,data=new_pca_frame,alternative='greater')


## Models of bacterial parameters, with random depth
## effect and linear model for model comparison using AIC
bac_modeling<-cbind(bulk_mat,bac_pca)
bac_gc_model_linear<-lm(data=bac_modeling,formula=gene_GC~nreads+fraction)
bac_nc_model_linear<-lm(data=bac_modeling,formula=N_C_rat~nreads+fraction)
bac_narsc_model_linear<-lm(data=bac_modeling,formula=N_ARSC~nreads+fraction)
bac_gc_model<-lme(fixed=gene_GC~fractions+log10(nreads),random=~1|depth,
                  data=bac_modeling,
                  cor=corExp(),
                  method='ML')
bac_nc_model<-lme(fixed=N_C_rat~fractions+log10(nreads),random=~1|depth,
                  data=bac_modeling,
                  cor=corExp(),
                  method='ML')
bac_narsc_model<-lme(fixed=N_ARSC~fractions+log10(nreads),random=~1|depth,
                     data=bac_modeling,
                     cor=corExp(),
                     method='ML')

## Archaeal parameters
arc_modeling<-cbind(bulk_mat,arc_pca)
arc_gc_model_linear<-lm(data=arc_modeling,formula=gene_GC~nreads+fraction)
arc_nc_model_linear<-lm(data=arc_modeling,formula=N_C_rat~nreads+fraction)
arc_narsc_model_linear<-lm(data=arc_modeling,formula=N_ARSC~nreads+fraction)
arc_gc_model<-lme(fixed=gene_GC~fractions+log10(nreads)+1,random=~1|depth,
                  data=arc_modeling,
                  cor=corExp(),method='ML')
arc_nc_model<-lme(fixed=N_C_rat~fractions+log10(nreads),random=~1|depth,
                  data=arc_modeling,
                  cor=corExp(),method='ML')
arc_narsc_model<-lme(fixed=N_ARSC~fractions+log10(nreads),random=~1|depth,
                     data=arc_modeling,
                     cor=corExp(),method='ML')

## Viral parameters

vir_modeling<-cbind(bulk_mat,viral_pca)
vir_gc_model_linear<-lm(data=vir_modeling,formula=gene_GC~nreads+fraction)
vir_nc_model_linear<-lm(data=vir_modeling,formula=N_C_rat~nreads+fraction)
vir_narsc_model_linear<-lm(data=vir_modeling,formula=N_ARSC~nreads+fraction)
vir_gc_model<-lme(fixed=gene_GC~fractions+log10(nreads),random=~1|depth,
                  data=vir_modeling,
                  cor=corExp(),method='ML')
vir_nc_model<-lme(fixed=N_C_rat~fractions+log10(nreads),random=~1|depth,
                  data=vir_modeling,
                  cor=corExp(),method='ML')
vir_narsc_model<-lme(fixed=N_ARSC~fractions+log10(nreads),random=~1|depth,
                     data=vir_modeling,
                     cor=corExp(),method='ML')

## Extracting the random effects
random_effect_frame<-data.frame(depth=as.numeric(rownames(bac_gc_model$coefficients$random$depth)),
                                bac_gc=bac_gc_model$coefficients$random$depth[,1],
                                bac_nc=bac_nc_model$coefficients$random$depth[,1],
                                bac_narsc=bac_narsc_model$coefficients$random$depth[,1],
                                arc_gc=arc_gc_model$coefficients$random$depth[,1],
                                arc_nc=arc_nc_model$coefficients$random$depth[,1],
                                arc_narsc=arc_narsc_model$coefficients$random$depth[,1],
                                vir_gc=vir_gc_model$coefficients$random$depth[,1],
                                vir_nc=vir_nc_model$coefficients$random$depth[,1],
                                vir_narsc=vir_narsc_model$coefficients$random$depth[,1])

## Plotting the random effects (supplementary figure)
re_bac_gc<-ggplot(random_effect_frame)+
  geom_point(aes(x=depth,
                 y=bac_gc))+
  geom_smooth(aes(x=depth,
                 y=bac_gc))+
  scale_x_reverse()+
  coord_flip()+
  xlab('Depth [m]')+
  ylab('Bacteria GC Random Effect')+
  theme_bw()+
  theme(text=element_text(size=14))
re_bac_nc<-ggplot(random_effect_frame)+
  geom_point(aes(x=depth,
                 y=bac_nc))+
  geom_smooth(aes(x=depth,
                 y=bac_nc))+
  scale_x_reverse()+
  coord_flip()+
  xlab('Depth [m]')+
  ylab('Bacteria NC Random Effect')+
  theme_bw()+
  theme(text=element_text(size=14))
re_bac_narsc<-ggplot(random_effect_frame)+
  geom_point(aes(x=depth,
                 y=bac_narsc))+
  geom_smooth(aes(x=depth,
                 y=bac_narsc))+
  scale_x_reverse()+
  coord_flip()+
  xlab('Depth [m]')+
  ylab('Bacteria N-ARSC Random Effect')+
  theme_bw()+
  theme(text=element_text(size=14))
re_arc_gc<-ggplot(random_effect_frame)+
  geom_point(aes(x=depth,
                 y=arc_gc))+
  geom_smooth(aes(x=depth,
                 y=arc_gc))+
  scale_x_reverse()+
  coord_flip()+
  xlab('Depth [m]')+
  ylab('Archaea GC Random Effect')+
  theme_bw()+
  theme(text=element_text(size=14))
re_arc_nc<-ggplot(random_effect_frame)+
  geom_point(aes(x=depth,
                 y=arc_nc))+
  geom_smooth(aes(x=depth,
                 y=arc_nc))+
  scale_x_reverse()+
  coord_flip()+
  xlab('Depth [m]')+
  ylab('Archaea NC Random Effect')+
  theme_bw()+
  theme(text=element_text(size=14))
re_arc_narsc<-ggplot(random_effect_frame)+
  geom_point(aes(x=depth,
                 y=arc_narsc))+
  geom_smooth(aes(x=depth,
                 y=arc_narsc))+
  scale_x_reverse()+
  coord_flip()+
  xlab('Depth [m]')+
  ylab('Archaea N-ARSC Random Effect')+
  theme_bw()+
  theme(text=element_text(size=14))
re_vir_gc<-ggplot(random_effect_frame)+
  geom_point(aes(x=depth,
                 y=vir_gc))+
  geom_smooth(aes(x=depth,
                 y=vir_gc))+
  scale_x_reverse()+
  coord_flip()+
  xlab('Depth [m]')+
  ylab('Virus GC Random Effect')+
  theme_bw()+
  theme(text=element_text(size=14))
re_vir_nc<-ggplot(random_effect_frame)+
  geom_point(aes(x=depth,
                 y=vir_nc))+
  geom_smooth(aes(x=depth,
                 y=vir_nc))+
  scale_x_reverse()+
  coord_flip()+
  xlab('Depth [m]')+
  ylab('Virus NC Random Effect')+
  theme_bw()+
  theme(text=element_text(size=14))
re_vir_narsc<-ggplot(random_effect_frame)+
  geom_point(aes(x=depth,
                 y=vir_narsc))+
  geom_smooth(aes(x=depth,
                 y=vir_narsc))+
  scale_x_reverse()+
  coord_flip()+
  xlab('Depth [m]')+
  ylab('Virus N-ARSC Random Effect')+
  theme_bw()+
  theme(text=element_text(size=14))

fig_s2<-(re_bac_gc|re_bac_nc|re_bac_narsc)/(plot_spacer()|re_arc_nc|re_arc_narsc)/(plot_spacer()|re_vir_nc|re_vir_narsc)
ggsave('figures/SF2.pdf',fig_s2)


## Producing model parameter summary supplemental table
full_table_list<-lapply(list(bac_gc_model,bac_nc_model,bac_narsc_model,
                             arc_gc_model,arc_nc_model,arc_narsc_model,
                             vir_gc_model,vir_nc_model,vir_narsc_model),
                        broom.mixed::tidy)
categories_1<-rep(c('Bacteria','Archaea','Virus'),each=3)
categories_2<-rep(c('Gene GC','Gene N:C Ratio','Gene #N per Side Chain'),3)
for(i in 1:length(full_table_list)){
  full_table_list[[i]]<-data.frame(Domain=categories_1[i],
                                   Parameter=categories_2[i],
                                   full_table_list[[i]])
}
full_table<-do.call(rbind,full_table_list) %>%
  rename('Effect'='effect',
         'Covariate'='term',
         'Estimate'='estimate',
         'Standard Error'='std.error',
         'Degrees of Freedom'='df',
         'Test Statistic'='statistic',
         'p'='p.value') %>%
  mutate(Covariate=factor(Covariate,
                          levels=c('(Intercept)',
                                   'fractionsSV',
                                   'log10(nreads)',
                                   'sd_(Intercept)',
                                   'sd_Observation'),
                          labels=c('Intercept','Planktonic Fraction',
                                   'Sequencing Depth (log10 reads)',
                                   'Depth Random Effect',
                                   'Intercept Random Effect')),
         Effect=ifelse(Effect=='fixed','Fixed','Random'),
         Covariate=as.character(Covariate),
         Estimate=signif(Estimate,digits=3),
         `Standard Error`=signif(`Standard Error`,digits=3),
         `Test Statistic`=signif(`Test Statistic`,digits=3),
         p=signif(p,digits=3)) %>%
  select(-group)

library(gt)
regression_table<-gt(full_table) %>%
  tab_header(title='Mixed Effect Model Summaries',
             subtitle='Test for Significant Effect of Size-Fraction on Stoichiogenomic Parameters') %>%
  tab_row_group(rows=31:45,
                label='Viral Models') %>%
  tab_row_group(rows=16:30,
                label='Archaeal Models') %>%
  tab_row_group(rows=1:15,
                label='Bacterial Models')

gtsave(regression_table,'figures/supplemental_table2.html')


