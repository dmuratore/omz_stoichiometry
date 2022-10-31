## Doing amino acid composition analysis on core genes for each
## Domain
library(glmnet)
library(tidyverse)
## Starting with viral Gp23
## Retrieve genes ids annotated as Gp23 from virus gene database
capsid_ids<-grep("Gp23",virus_information[[1]]$description)
## Retrieve the converages of those genes from all samples
capsid_cov<-virus_information[[2]][capsid_ids,]
## Get the amino acid abundances for each gene
capsid_traits<-virus_information[[1]][capsid_ids,c(2,3,5:29)]
## Construct weighted average of the coverages of each gene for each amino acid
relative_capsid<-colSums(t(t(capsid_cov)/virus_information[[3]]))
## Normalize
capsid_rel_abund<-apply(capsid_cov,2,function(x) x/sum(x))
weighted_capsid_traits<-data.frame(t(capsid_rel_abund)%*%as.matrix(capsid_traits)) %>%
  rownames_to_column(var='sample')
## Join with contextual data
capsid_frame<-full_join(mutate(bulk_mat,sample=paste0(sample,'_cov')),weighted_capsid_traits) %>%
  full_join(rownames_to_column(data.frame(relative_capsid),var='sample'))
## Rescale data for LASSO implementation
capsid_a<-capsid_frame %>% select(c(N_ARSC,N_C_rat,fraction)) %>% 
  drop_na() %>%
  mutate(fraction=ifelse(fraction=='PF',1,0)) %>%
  mutate(N_ARSC=scale(N_ARSC),
         N_C_rat=scale(N_C_rat))
capsid_b<-scale(capsid_frame[,19:38])
## Split into test train and learn model
set.seed(175)
v_train<-sample(1:nrow(capsid_a),size=round(0.7*nrow(capsid_a)))
set.seed(175)
virus_cv<-cv.glmnet(as.matrix(capsid_b[v_train,]),
                    as.matrix(capsid_a[v_train,]),
                    family='mgaussian')
## Pull out optimal model using CV-determined L1-norm penalty
set.seed(175)
virus_lasso<-glmnet(capsid_b[v_train,],
                    capsid_a[v_train,],
                    family='mgaussian',
                    lambda=virus_cv$lambda.1se)
## Retrieve model parameters and assess model fit
v_pars<-cbind(rownames(as.matrix(virus_lasso$beta$N_C_rat)),
              as.matrix(virus_lasso$beta$N_C_rat))
virus_test<-assess.glmnet(virus_lasso,
                          newx=as.matrix(capsid_b[-v_train,]),
                          newy=as.matrix(capsid_a[-v_train,]))

## Move on to archaea, same procedure
arc_marker<-grep("K03531",arc_information[[1]]$kegg_orths)
arc_cov<-arc_information[[2]][arc_marker,]
arc_traits<-arc_information[[1]][arc_marker,c(2,3,5:29)]
arc_rel_abund<-apply(arc_cov,2,function(x) x/sum(x))
weighted_arc_traits<-data.frame(t(arc_rel_abund)%*%as.matrix(arc_traits)) %>%
  rownames_to_column(var='sample')
arc_frame<-full_join(mutate(bulk_mat,sample=paste0(sample,'_cov')),weighted_arc_traits)
arc_a<-arc_frame %>% select(c(N_ARSC,N_C_rat,fraction)) %>% 
  drop_na() %>%
  mutate(fraction=ifelse(fraction=='PF',1,0)) %>%
  mutate(N_ARSC=scale(N_ARSC),
         N_C_rat=scale(N_C_rat))
arc_b<-scale(arc_frame[,19:38] %>% drop_na())
set.seed(25739)
a_train<-sample(1:nrow(arc_a),size=round(0.7*nrow(arc_a)))
set.seed(25739)
arc_cv<-cv.glmnet(as.matrix(arc_b[a_train,]),
                  as.matrix(arc_a[a_train,]),
                  family='mgaussian')
set.seed(25739)
arc_lasso<-glmnet(arc_b[a_train,],
                  arc_a[a_train,],
                  family='mgaussian',
                  lambda=arc_cv$lambda.1se)
a_pars<-cbind(rownames(as.matrix(arc_lasso$beta$N_C_rat)),
              as.matrix(arc_lasso$beta$N_C_rat))
cbind(arc_lasso$beta$fraction@Dimnames[[1]][arc_lasso$beta$fraction@i[-1]],
      arc_lasso$beta$fraction@x[-1])
arc_test<-assess.glmnet(arc_lasso,
                        newx=as.matrix(arc_b[-a_train,]),
                        newy=as.matrix(arc_a[-a_train,]))

## Finally bacteria
bac_marker<-grep("K03060",bac_information[[1]]$kegg_orths)
bac_cov<-bac_information[[2]][bac_marker,]
bac_traits<-bac_information[[1]][bac_marker,c(2,3,5:29)]
bac_rel_abund<-apply(bac_cov,2,function(x) x/sum(x))
weighted_bac_traits<-data.frame(t(bac_rel_abund)%*%as.matrix(bac_traits)) %>%
  rownames_to_column(var='sample')
bac_frame<-full_join(mutate(bulk_mat,sample=paste0(sample,'_cov')),weighted_bac_traits)
bac_a<-bac_frame %>% select(c(N_ARSC,N_C_rat,fraction)) %>% 
  drop_na() %>%
  mutate(fraction=ifelse(fraction=='PF',1,0)) %>%
  mutate(N_ARSC=scale(N_ARSC),
         N_C_rat=scale(N_C_rat))
bac_b<-scale(bac_frame[,19:38] %>% drop_na())
set.seed(1109)
b_train<-sample(1:nrow(bac_a),size=round(0.7*nrow(bac_a)))
set.seed(1109)
bac_cv<-cv.glmnet(as.matrix(bac_b[b_train,]),
                  as.matrix(bac_a[b_train,]),
                  family='mgaussian')
set.seed(1109)
bac_lasso<-glmnet(bac_b[b_train,],
                  bac_a[b_train,],
                  family='mgaussian',
                  lambda=bac_cv$lambda.1se)
b_pars<-cbind(rownames(as.matrix(bac_lasso$beta$N_C_rat)),
              as.matrix(bac_lasso$beta$N_C_rat))
bac_test<-assess.glmnet(bac_lasso,
                        newx=as.matrix(bac_b[-b_train,]),
                        newy=as.matrix(bac_a[-b_train,]))

## Visualizing model performance on test data
par(mfrow=c(3,1))
pdf(file='figures/SF4.pdf')
plot(predict(arc_lasso,as.matrix(arc_b[-a_train,]))[,2,1],
     arc_a[-a_train,2],
     main='Archaeal ftsZ Amino Acid Model Test\nSlope=1.0072, R^2=99.42%',
     ylab='Actual N:C Ratio',
     xlab='Predicted N:C Ratio')
abline(b=1,a=0)
summary(lm(arc_a[-a_train,2]~predict(arc_lasso,as.matrix(arc_b[-a_train,]))[,2,1]))

plot(predict(bac_lasso,as.matrix(bac_b[-b_train,]))[,2,1],
     bac_a[-b_train,2],
     main='Bacterial rpoZ Amino Acid Model Test\nSlope=0.961, R^2=93.24%',
     ylab='Actual N:C Ratio',
     xlab='Predicted N:C Ratio')
abline(b=1,a=0)
summary(lm(bac_a[-b_train,2]~predict(bac_lasso,as.matrix(bac_b[-b_train,]))[,2,1]))

plot(predict(virus_lasso,as.matrix(capsid_b[-v_train,]))[,2,1],
     capsid_a[-v_train,2],
     main='Virus Gp23 Amino Acid Model Test\nSlope=0.918, R^2=92.22%',
     ylab='Actual N:C Ratio',
     xlab='Predicted N:C Ratio')
abline(b=1,a=0)
summary(lm(capsid_a[-v_train,2]~predict(virus_lasso,as.matrix(capsid_b[-v_train,]))[,2,1]))
dev.off()
## Generating figure 4 from model parameters
model_pars<-rbind(b_pars,a_pars,v_pars)
model_pars<-data.frame(aa=model_pars[,1],
                       par=as.numeric(model_pars[,2]),
                       org=c(rep('Bacteria',nrow(b_pars)),
                             rep('Archaea',nrow(a_pars)),
                             rep('Virus',nrow(v_pars)))) %>%
  mutate(name=factor(aa,levels=c('C','E','D','G',
                                 'H','K','L','N',
                                 'Q','P','T','I',
                                 'S','R','V','A',
                                 'F','M','W','Y'),
                     labels=c('Cys','Glu',
                              'Asp','Gly',
                              'His','Lys',
                              'Leu','Asn',
                              'Gln','Pro',
                              'Thr','Ile',
                              'Ser','Arg',
                              'Val','Ala',
                              'Phe','Met','Trp','Tyr')),
         n_arsc=ifelse(aa=='R',3,
                       ifelse(aa=='H',2,
                              ifelse(aa %in% c('K','N','Q','W'),
                                     1,0)))) %>%
  mutate(name=tidytext::reorder_within(name,par,org)) %>%
  filter(par!=0)
lasso_model<-ggplot(model_pars)+
  geom_segment(aes(x=0,y=name,
                   yend=name,
                   xend=par,
                   col=factor(n_arsc)),
               size=2)+
  facet_wrap(~factor(org,
                     levels=c('Bacteria','Archaea','Virus'),
                     labels=c('Bacterial rpoZ','Archaeal ftsZ','Viral Gp23')),scales='free',ncol=1)+
  tidytext::scale_y_reordered()+
  xlab('Coefficient of Amino Acid on Protein N:C Ratio')+
  scale_color_viridis_d(name='# N Atoms in Side Chain')+
  ylab('')+
  theme_bw()+
  theme(strip.background=element_blank(),
        strip.text=element_text(size=12))
ggsave('figures/F4.pdf',device='pdf',scale=1.25)

## Writing model output to table
gt::gt(model_pars %>% group_by(org) %>% arrange(desc(par),desc(n_arsc))) %>%
  tab_header(
    title='Amino Acid Composition Models for Core Gene N:C Ratios'
  ) %>%
  cols_move_to_start(c(name,aa)) %>%
  cols_label(name='Amino Acid/Domain',
             aa='Amino Acid',
             par='N:C Effect',
             n_arsc='N-ARSC of Amino Acid') %>%
  gtsave('figures/supplemental_table3.html')

## Generating models for depth and size-fraction distributions
## for core genes

## Model with depth as random effect
bac_marker_model<-lme(fixed=N_C_rat~fraction+log10(nreads),random=~1|depth,
                      data=drop_na(bac_frame),
                      cor=corExp(),
                      method='ML')
vir_marker_model<-lme(fixed=N_C_rat~log10(nreads)+fraction,random=~1|depth,
                      data=capsid_frame,
                      cor=corExp(),
                      method='ML')
arc_marker_model<-lme(fixed=N_C_rat~fraction+log10(nreads),random=~1|depth,
                      data=drop_na(arc_frame),
                      cor=corExp(),
                      method='ML')




## Model without depth effect for model comparison/likelihood ratio testing
bac_marker_linear<-lm(data=drop_na(bac_frame),formula=N_C_rat~nreads+fraction)
vir_marker_linear<-lm(data=capsid_frame,formula=N_C_rat~log10(nreads)+fraction)
arc_marker_linear<-lm(data=drop_na(arc_frame),formula=N_C_rat~log10(nreads)+fraction)

## Modeling the correlation of depth effect with depth
cor.test(as.numeric(rownames(bac_marker_model$coefficients$random$depth)),
         bac_marker_model$coefficients$random$depth[,1],method='spearman')
cor.test(as.numeric(rownames(arc_marker_model$coefficients$random$depth)),
         arc_marker_model$coefficients$random$depth[,1],method='spearman')
cor.test(as.numeric(rownames(vir_marker_model$coefficients$random$depth)),
         vir_marker_model$coefficients$random$depth[,1],method='spearman')

## Model comparisons 
anova(bac_marker_model,bac_marker_linear)
anova(vir_marker_model,vir_marker_linear)
anova(arc_marker_model,arc_marker_linear)

## Producing model parameter summary supplemental table
full_table_list<-lapply(list(bac_marker_model,
                             vir_marker_model,
                             arc_marker_model),
                        broom.mixed::tidy)
categories_1<-rep(c('Bacteria','Archaea','Virus'),each=1)
categories_2<-rep(c('rpoZ N:C','ftsZ N:C','Gp23 N:C'),1)
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
                                   'fractionSV',
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
             subtitle='Test for Significant Effect of Size-Fraction on Marker Genes') %>%
  tab_row_group(rows=11:15,
                label='Viral Models') %>%
  tab_row_group(rows=6:10,
                label='Archaeal Models') %>%
  tab_row_group(rows=1:5,
                label='Bacterial Models')

gtsave(regression_table,'figures/supplemental_table4.html')

## Plotting data
mini_bac_frame<-pivot_longer(bac_frame,names_to='amino_acid',
                             cols=c('W','C','H','M','K','L',
                                    'A','E','I','F','Y','V',
                                    'D','R','S','T','N','Q','G','P')) %>%
  mutate(aa_polarity=ifelse(amino_acid %in% c('A','I','L','V','F','W','C','M'),
                            'hydrophobic',
                            ifelse(amino_acid %in% c('P','H','Y','T','S'),
                                   'neutral','hydrophilic')))

b_nc<-ggplot(mini_bac_frame)+
  geom_point(aes(x=depth,y=N_C_rat,col=fraction))+
  scale_x_reverse(lim=c(300,0))+
  scale_y_continuous()+
  coord_flip()+
  geom_smooth(aes(x=depth,y=N_C_rat,col=fraction))+
  theme_bw()+
  theme(text=element_text(size=14))+
  scale_color_manual(name='Size Fraction',
                     values=c('darkgreen','cyan3'),
                     labels=c('Particulate','Free-Living'))+
  ylab('Bacterial rpoZ Protein Sequence N:C Ratio')+
  xlab('Depth [m]')

mini_arc_frame<-pivot_longer(arc_frame,names_to='amino_acid',
                             cols=c('W','C','H','M','K','L',
                                    'A','E','I','F','Y','V',
                                    'D','R','S','T','N','Q','G','P')) %>%
  mutate(aa_polarity=ifelse(amino_acid %in% c('A','I','L','V','F','W','C','M'),
                            'hydrophobic',
                            ifelse(amino_acid %in% c('P','H','Y','T','S'),
                                   'neutral','hydrophilic')))
a_nc<-ggplot(mini_arc_frame)+
  geom_point(aes(x=depth,y=N_C_rat,col=fraction))+
  scale_x_reverse(lim=c(300,0))+
  scale_y_continuous()+
  coord_flip()+
  geom_smooth(aes(x=depth,y=N_C_rat,col=fraction))+
  scale_color_manual(name='Size Fraction',
                     values=c('darkgreen','cyan3'),
                     labels=c('Particulate','Free-Living'))+
  theme_bw()+
  theme(text=element_text(size=14))+
  ylab('Archaeal ftsZ Protein Sequence N:C Ratio')+
  xlab('Depth [m]')

mini_capsid_frame<-pivot_longer(capsid_frame,names_to='amino_acid',
                                cols=c('W','C','H','M','K','L',
                                       'A','E','I','F','Y','V',
                                       'D','R','S','T','N','Q','G','P')) %>%
  mutate(aa_polarity=ifelse(amino_acid %in% c('A','I','L','V','F','W','C','M'),
                            'hydrophobic',
                            ifelse(amino_acid %in% c('P','H','Y','T','S'),
                                   'neutral','hydrophilic')))

v_nc<-ggplot(mini_capsid_frame)+
  geom_point(aes(x=depth,y=N_C_rat))+
  scale_x_reverse(lim=c(300,0))+
  scale_y_continuous(lim=c(0.1,0.125))+
  coord_flip()+
  geom_smooth(aes(x=depth,y=N_C_rat),col='gold4')+
  theme_bw()+
  theme(text=element_text(size=14))+
  ylab('Capsid Protein Sequence N:C Ratio')+
  xlab('Depth [m]')

supp_figure<-b_nc|a_nc|v_nc
ggsave(supp_figure,filename='figures/SF3.pdf',scale=2)

