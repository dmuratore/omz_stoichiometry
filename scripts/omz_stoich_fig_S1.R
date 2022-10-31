## Taxonomic Profiling from Metaxa2 reads
## All necessary files and data to reproduce this analysis are
## zipped together in the 16s zip in the data subdirectory
library(dada2)
library(msa)
library(phyloseq)
library(data.table)
library(tidyverse)
library(RColorBrewer)
setwd('data/16s')

## Get all forward and reverse read file names
forward_reads<-sort(list.files('f_reads',full.names=TRUE))
reverse_reads<-sort(list.files('r_reads',full.names=TRUE))

filter_forward_outputnames<-file.path('filtered_reads/forward',paste0(gsub('.fq','',list.files('f_reads/')),'_filtered.fq'))
filter_reverse_outputnames<-file.path('filtered_reads/reverse',paste0(gsub('.fq','',list.files('r_reads/')),'_filtered.fq'))

## Filter and trimming
out<-filterAndTrim(forward_reads,
                   filter_forward_outputnames,
                   rev=reverse_reads,
                   filt.rev=filter_reverse_outputnames,
                   compress=TRUE,
                   multithread=TRUE,maxEE=c(2,2),
                   truncQ=2)

## Learning per base transition error rates
f_errs<-learnErrors(filter_forward_outputnames,
                    multithread=TRUE,
                    verbose=1,
                    randomize=TRUE)
r_errs<-learnErrors(filter_reverse_outputnames,
                    multithread=TRUE,
                    verbose=1,
                    randomize=TRUE)

## Dereplicate and profile
derep_f<-derepFastq(filter_forward_outputnames,n=1e8,verbose=TRUE)
derep_r<-derepFastq(filter_reverse_outputnames,n=1e8,verbose=TRUE)
dada_f<-dada(derep_f,f_errs,multithread=TRUE,verbose=1)
dada_r<-dada(derep_r,r_errs,multithread=TRUE,verbose=1)

## Merging and tabularizing results
merged<-mergePairs(dada_f,derep_f,dada_r,derep_r,
                   minOverlap=6,verbose=TRUE)
seq_tab<-makeSequenceTable(merged)

## Assign taxonomy using silva database
taxa<-assignTaxonomy(seq_tab,
                     refFasta='silva_nr_v132_train_set.fa.gz',
                     multithread = TRUE,
                     tryRC=TRUE)

## Reorganizing taxonomy to split up proteobacteria
kingdom_annotes<-sapply(unique(taxa[,1]),function(x) length(which(taxa[,1]==x)))
phy_annotes<-sapply(unique(taxa[,2]),function(x) length(which(taxa[,2]==x)))

tax_with_class<-taxa
tax_with_class[which(tax_with_class[,2]=='Proteobacteria'),2]<-tax_with_class[which(tax_with_class[,2]=='Proteobacteria'),3]

## Getting environmental covariate data
sample<-rownames(seq_tab)
year<-gsub('ETNP','',sample)
year<-gsub('S[0-9].*$','',year)
stn<-gsub('ETNP[0-9][0-9]','',sample)
stn<-gsub('M.*$','',stn)
stn<-gsub('D.*$','',stn)
depth<-gsub('^.*G','',sample)
depth<-as.numeric(gsub('[A-Z].*$','',depth))
fraction<-gsub('^.*G[0-9]*','',sample)
fraction<-gsub('\\..*$','',fraction)
sample_frame<-data.frame(sample,year,stn,depth,fraction)
rownames(sample_frame)<-sample

## Incorporate taxonomy and sequence count table
phylo_object<-phyloseq(otu_table(seq_tab,taxa_are_rows = FALSE),
                       tax_table(tax_with_class),
                       sample_data(sample_frame))

## Summarize to Phylum level
glommed_to_phylum<-tax_glom(phylo_object,taxrank='Phylum')

## Data reformatting for plotting
phylum_counts<-glommed_to_phylum@otu_table@.Data
phylum_sample_data<-as.data.frame(glommed_to_phylum@sam_data@.Data)
colnames(phylum_sample_data)<-c('sample','year','station','depth','fraction')
phylum_sample_data$station<-gsub('MG.*$','',gsub('ETNP1[0-9]','',phylum_sample_data$sample))
phylum_phy_names<-glommed_to_phylum@tax_table@.Data
colnames(phylum_counts)<-phylum_phy_names[match(colnames(phylum_counts),rownames(phylum_phy_names)),2]

## Phyla to emphasize for plotting
keeper_phy<-c('Alphaproteobacteria',
              'Deltaproteobacteria',
              'Dinoflagellata','Protalveolata',
              'Thaumarchaeota','Cyanobacteria',
              'Marinimicrobia_(SAR406_clade)',
              'Gammaproteobacteria','Euryarchaeota',
              'Actinobacteria','Bacteroidetes','Verrucomicrobia',
              'Nanoarchaeaeota','Nitrospinae',
              'Epsilonbacteraeota','Acidobacteria')

## Creating bar chart counting frame
phylum_16s_counting_frame<-cbind(phylum_counts,phylum_sample_data) %>%
  gather(key='phylum',value='count',colnames(phylum_counts)) %>%
  group_by(sample) %>%
  mutate(rel_count=count/sum(count)) %>%
  mutate(kingdom=phylum_phy_names[match(phylum,phylum_phy_names[,2]),1]) %>%
  mutate(depth_num=gsub('_.*$','',gsub('^.*MG','',sample))) %>%
  group_by(phylum) %>%
  mutate(update_phy=ifelse(phylum %in% keeper_phy,phylum,'Other'))
phylum_16s_counting_frame$depth_num<-factor(phylum_16s_counting_frame$depth_num,
                                            levels=unique(phylum_16s_counting_frame$depth_num[order(phylum_16s_counting_frame$depth,decreasing=TRUE)]))

## Making supplemental plot
figure_s1<-phylum_16s_counting_frame %>%
  ggplot(aes(x=depth_num,y=count,fill=update_phy),col='black')+
  geom_col()+
  facet_wrap(~station*fraction,scales='free')+
  theme(legend.position='bottom',
        strip.background=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(angle=60,vjust=0.5))+
  coord_flip()+
  scale_x_discrete(name='Depth')+
  scale_y_continuous(name='Number 16S Reads')+
  scale_fill_manual(name='Phylum',
                    guide=guide_legend(nrow=4),
                    values=colorRampPalette(brewer.pal(12,"Set3"))(length(unique(phylum_16s_counting_frame$update_phy))))+
  theme(text=element_text(size=18))
ggsave(filename='../../figures/SF1.pdf',figure_s1,
       device='pdf',units='in',height=12,width=12)
