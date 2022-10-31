## Analyzing Virus Amino Acid Frequencies
## Reformatting the data used in the figure 4 script
mini_capsid_frame<-pivot_longer(capsid_frame,names_to='amino_acid',
                                cols=c('W','C','H','M','K','L',
                                       'A','E','I','F','Y','V',
                                       'D','R','S','T','N','Q','G','P')) %>%
  mutate(aa_polarity=ifelse(amino_acid %in% c('A','I','L','V','F','W','C','M'),
                            'hydrophobic',
                            ifelse(amino_acid %in% c('P','H','Y','T','S'),
                                   'neutral','hydrophilic')))
## Making plot
aa_traces<-ggplot(filter(mini_capsid_frame,amino_acid %in% c('R','Q','H','N','M','C')) %>%
         mutate(amino_acid=factor(amino_acid,levels=c('R','Q','N','H','M','C'))))+
  geom_point(aes(x=depth,y=value),size=3)+
  geom_smooth(aes(y=value,x=depth,group=amino_acid,col=amino_acid),
              method='gam')+
  facet_wrap(vars(amino_acid,aa_polarity),scale='free_x',nrow=2,
             labeller=as_labeller(c('W'='Tryptophan',
                                    'C'='Cysteine',
                                    'H'='Histidine',
                                    'M'='Methionine',
                                    'K'='Lysine',
                                    'L'='Leucine',
                                    'A'='Alanine',
                                    'E'='Glutamic Acid',
                                    'I'='Isoleucine',
                                    'F'='Phenylalanine',
                                    'Y'='Tyrosine',
                                    'V'='Valine',
                                    'D'='Aspartic Acid',
                                    'R'='Arginine',
                                    'S'='Serine',
                                    'T'='Threonine',
                                    'N'='Asparagine',
                                    'Q'='Glutamine',
                                    'G'='Glycine',
                                    'P'='Proline',
                                    'neutral'='Neutral',
                                    'hydrophobic'='Hydrophobic',
                                    'hydrophilic'='Polar')))+
  scale_x_reverse(lim=c(500,0))+
  theme_bw()+
  xlab('Depth [m]')+
  ylab('Avg # AA per\nCapsid Protein Sequence')+
  scale_color_brewer(palette='Set2',name='Amino Acid')+
  coord_flip()+
  theme(text=element_text(size=16),
        strip.background=element_rect(fill='white'),
        legend.position='none')

## Saving plot
ggsave('figures/F5.pdf',aa_traces,device='pdf',scale=1.25)

## Calculating AA frequency correlations
virus_ctd<-left_join(capsid_frame %>% select(-names(v_pars[which(as.numeric(v_pars[,2])<0),1])),nutrient_and_ctd %>% select(-year),
                     by=c('stn'='updated_station',
                          'depth'='Depth..m.'))
capsid_aa_correlations<-apply(virus_ctd[,19:29],2,function(x) cor.test(x,virus_ctd$oxygen,
                                                                       method='spearman'))
capsid_aa_stats<-do.call(rbind,lapply(capsid_aa_correlations,function(x) c(x$estimate,
                                                                           x$statistic,
                                                                           x$p.value))) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename('test_statistic'='S','p_value'='V3')
