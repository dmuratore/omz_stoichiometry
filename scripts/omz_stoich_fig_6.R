## Functional gene analysis
library(tidyverse)
library(broom)
## Using sdhC as our baseline gene for log-ratio differential abundance
## analyses
baseline_ko<-'ko:K00241' #sdhC

## Extracting pathway and annotation information
pathway_identifying<-function(row){
  kos<-strsplit(row[2],',')[[1]]
  no_map<-kos[grep('ko',kos)]
  if(length(kos)>0){
    output<-cbind(rep(row[1],length(no_map)),no_map)
  }else return(NULL)
  return(output)
}

bac_gene_to_out<-cbind(bac_information[[1]]$gene.id,bac_information[[1]]$kegg_orths)

## Reformatting data for analysis
gene_pathway_table<-apply(bac_gene_to_out,1,pathway_identifying)
n_paths<-do.call(c,lapply(gene_pathway_table,nrow))
gene_pathway_table<-do.call(rbind,gene_pathway_table)

## Preparing data structures to fill in with coverage data of KOs
coverage_pathway_table<-list()

for(i in 1:length(n_paths)){
  coverage_pathway_table[[i]]<-matrix(bac_information[[2]][i,],
                                      ncol=ncol(bac_information[[2]]),
                                      nrow=n_paths[i],byrow=TRUE)
}

coverage_pathway_table<-do.call(rbind,coverage_pathway_table)

## Uniting coverage and KO information
coverage_frame<-data.frame(gene_pathway_table,coverage_pathway_table)
## Normalizing coverages
coverage_table<-coverage_frame %>%
  split(coverage_frame$no_map) %>%
  lapply(.,function(x) colSums(matrix(unlist(x[,-c(1,2)]),ncol=58))) %>%
  do.call(rbind,.)

## Converting to log ratios and reformatting data for modeling
ratio_mat<-t(apply(coverage_table,1,function(x) x/coverage_table[baseline_ko,]))
sample_names<-colnames(bac_information[[2]])
colnames(ratio_mat)<-sample_names

## Incorporating sample covariates
sample_covars<-bac_information[[6]][,c('year','fractions','deeps')] %>%
  rownames_to_column(var='sample') %>%
  mutate(station=gsub('ETNP[0-9][0-9]','',gsub('MG.*$','',sample)),
         deeps=as.numeric(as.character(deeps)),
         year=as.character(year),
         station=gsub('f','',gsub('D.*$','',station))) %>%
  left_join(nutrient_and_ctd,by=c('station'='updated_station','year'='year','deeps'='Depth..m.'))

## Reformatting sample covariates for easier modeling syntax
full_ratio_mat<-data.frame(t(rbind(t(sample_covars),ratio_mat))) %>%
  gather('ko','ratio',contains('ko')) %>%
  mutate(ratio=as.numeric(ratio),
         deeps=as.numeric(as.character(deeps)),
         round_deep=round(deeps,-1),
         zone=if_else(deeps>80,if_else(deeps>500,'below','omz'),'above'))

## Removing samples wiht no coverage
full_ratio_mat_nozero<-full_ratio_mat %>%
  filter(!(ko %in% full_ratio_mat$ko[which(full_ratio_mat$ratio==0)])) %>%
  filter(ko!=gsub(':','.',baseline_ko)) %>%
  mutate(lograt=log(ratio))

## Applying linear models
model_list<-full_ratio_mat_nozero %>%
  split(full_ratio_mat_nozero$ko) %>%
  lapply(.,function(x) broom::tidy(lm(lograt~fractions*round_deep,data=x)))

## Extracting model parameters
pars<-do.call(rbind,lapply(model_list,function(x) x[c(2,3,4),c(1,2,4,5)])) %>%
  rownames_to_column() %>%
  arrange(p.value)

## Multiple testing control sequence
# Set FDR level let's say 0.10
q<-0.1
# Set number tests
m<-nrow(pars)
# Set ids
i<-1:m
# Set min thresholds
min_thresh<-i*q/m
# Find which are less
pass_thresh<-which(pars$p.value<=min_thresh)
# 87 Parameters pass initial threshold
# Calculate the slopes
slopes<-(1-pars$p.value)/(m+1-i)
first_fail<-min(which(diff(slopes)<0))
# Set m0hat
m0hat<-min((1/slopes[first_fail]+1),m)
# Compare p values
passes_m0hat<-min(which(!pars$p.value<=(i*q/m0hat)))
# Reject all k hypotheses
rejects<-pars[1:passes_m0hat,] %>%
  mutate(ko=gsub('ko','K',gsub('\\.[0-9]$','',rowname)))
length(unique(rejects$ko))

rejects %>% group_by(term) %>% count()
filter(rejects,term=='fractionsSV')$ko

## Pulling KO annotations from KEGG
sig_kos<-rejects$ko

short_names<-c('secF',
               'hpaI',
               'guaB',
               'kdsB',
               'livM',
               'potC',
               'rhmA',
               'HMGCS',
               'flgH',
               'uncharacterized',
               'pepN',
               'secD',
               'argF',
               'rpsJ',
               'MGAT3',
               'csrA',
               'yydH',
               'secDF',
               'gyaR',
               'atzF',
               'flgG',
               'smc',
               'lptD',
               'ilvI',
               'gcvP',
               'argB',
               'spoT',
               'nadE',
               'gcvPA',
               'kipI',
               'ridA',
               'tktA',
               'glpK',
               'uncharacterized',
               'pxpA',
               'pyrR',
               'fuvX',
               'thiG',
               'bkdA',
               'sdhB',
               'livG',
               'norM',
               'MVK',
               'pilB',
               'miaA',
               'plsX',
               'nadE',
               'garL',
               'TSFM',
               'bchG',
               'GALT29A',
               'protein-tyrosine phosphatase',
               'tig',
               'ycjG',
               'infB',
               'ALDH',
               'mlaA',
               'uncharacterized',
               'uncharacterized',
               'mmsA',
               'dacA',
               'trkA',
               'inlA',
               'modA',
               'pepP',
               'malZ',
               'gpmA',
               'prtC',
               'MGMT')

## Organizing into groups based on functions from KEGG pahtways
groups<-c('Protein export',
          'Amino acid synthesis',
          'Nucleotide metabolism',
          'Nucleotide metabolism',
          'Amino acid uptake',
          'Amino acid uptake',
          'Carbon metabolism',
          'Amino acid degradation',
          'Flagellar',
          'Other',
          'Amino acid degradation',
          'Protein export',
          'Arginine biosynthesis',
          'Other',
          'Carbon metabolism',
          'Carbon metabolism',
          'Amino acid degradation',
          'Protein export',
          'Other',
          'Arginine biosynthesis',
          'Flagellar',
          'Other',
          'Other',
          'Amino acid synthesis',
          'Amino acid degradation',
          'Arginine biosynthesis',
          'Nucleotide metabolism',
          'Other',
          'Amino acid degradation',
          'Other',
          'Other',
          'Amino acid synthesis',
          'Other',
          'Other',
          'Other',
          'Nucleotide metabolism',
          'Other',
          'Other',
          'Amino acid degradation',
          'Carbon metabolism',
          'Amino acid uptake',
          'Other',
          'Other',
          'Flagellar',
          'Other',
          'Other',
          'Other',
          'Carbon metabolism',
          'Other',
          'Other',
          'Other',
          'Amino acid synthesis',
          'Other',
          'Amino acid degradation',
          'Other',
          'Amino acid synthesis',
          'Other',
          'Other',
          'Other',
          'Amino acid degradation',
          'Other',
          'Other',
          'Other',
          'Other',
          'Amino acid degradation',
          'Carbon metabolism',
          'Carbon metabolism',
          'Other','Carbon metabolism')

## Preparing for plotting
rejection_rows<-coverage_table[gsub('\\.',':',gsub('\\.[0-9]$','',rejects$rowname)),]
colnames(rejection_rows)<-colnames(ratio_mat)

coverage_plotting_df<-cbind(rejects,
                            rejection_rows,short_names,groups) %>% 
  gather('sample','coverage',contains('ETNP')) %>%
  mutate(fraction=gsub('[0-9]','',gsub('\\.[0-9].*cov$','',gsub('^.*MG','',sample))),
         depth=as.numeric(gsub('[A-Z].*$','',gsub('^.*MG','',sample))),
         zone=if_else(depth>80,if_else(depth>500,'below','omz'),'above'),
         ko=gsub('^K','ko',ko)) %>%
  left_join(full_ratio_mat_nozero,by=c('sample'='sample','ko'='ko'))

## Creating figure
functional_analysis<-ggplot(coverage_plotting_df)+
  geom_jitter(aes(y=lograt,
                  x=fct_reorder(short_names,
                                desc(estimate)),
                  col=fractions),width=0.2,alpha=0.5)+
  coord_flip()+
  facet_wrap(~factor(groups,
                     levels=c('Nucleotide metabolism',
                              'Amino acid uptake','Amino acid degradation',
                              'Amino acid synthesis','Arginine biosynthesis',
                              'Protein export','Carbon metabolism',
                              'Flagellar','Other')),
             ncol=3,
             scale='free_y')+
  scale_color_manual(name='Size Fraction',
                     values=c('darkgreen','cyan3'),
                     labels=c('Particulate','Free-Living'))+
  theme_bw()+
  ylab('Gene Abundance (log-ratio)')+
  theme(axis.title.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=14))

ggsave('figures/F6.pdf',
       functional_analysis,scale=1.25)



