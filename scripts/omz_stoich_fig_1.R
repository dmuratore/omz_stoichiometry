## Code for creating figure 1 and conducting CTD data processing
## Metagenome-wide GC content and other statistics calculated using
## python code file gc_content.py

## Loading libraries
# sjPlot library from CRAN also required to generate table S1 html output
library(tidyverse)
library(reticulate)
library(nlme)
use_condaenv('base')
gsw<-import('gsw')

## Functions
ctd_process_2<-function(cast){
  if('Temperature.ITS90...C' %in% colnames(cast)){
    if('Fluorescence..ug.l.' %in% colnames(cast)){
      colnames(cast)[which(colnames(cast)=='Fluorescence..ug.l.')]<-'Fluorescence'
    }
    out<-cast %>%
      mutate(corrected_depth=gsw$z_from_p(`Depth..dB`,Latitude),
             sa=gsw$SA_from_SP(`Salinity..psu`,`Depth..dB`,Longitude,Latitude),
             ct=gsw$CT_from_t(sa,`Temperature.ITS90...C`,`Depth..dB`),
             sigma_theta=gsw$sigma0(sa,ct),
             o2_sol=gsw$O2sol(sa,ct,`Depth..dB`,Longitude,Latitude),
             o2_sat=100*(`Oxygen.SBE43..µmol.kg.1`/o2_sol),
             rounded_depth=-1*round((corrected_depth/1))*1) %>%
      rename(lat=Latitude,
             long=Longitude,
             temp=`Temperature.ITS90...C`,
             oxygen=`Oxygen.SBE43..µmol.kg.1`,
             fl=Fluorescence,
             pressure=`Depth..dB`,
             salinity=`Salinity..psu`) %>%
      select(lat,long,temp,
             salinity,sigma_theta,sa,
             ct,corrected_depth,pressure,
             oxygen,fl,o2_sol,o2_sat,rounded_depth) %>%
      group_by(rounded_depth) %>%
      summarize(across(everything(),mean,.names='{.col}')) %>%
      mutate(mld=rounded_depth[min(which(sigma_theta>=(sigma_theta[which(rounded_depth==10)]+0.125)))])
  }
  if('Oxygen.uM' %in% colnames(cast)){
    out<-cast %>%
      mutate(corrected_depth=gsw$z_from_p(`Pressure..db.`,Lat),
             sa=gsw$SA_from_SP(`Salinity`,`Pressure..db.`,Long,Lat),
             ct=gsw$CT_from_t(sa,`Temp.oC`,`Pressure..db.`),
             sigma_theta=gsw$sigma0(sa,ct),
             o2_sol=gsw$O2sol(sa,ct,`Pressure..db.`,Long,Lat),
             o2_sat=100*(`Oxygen.uM`/o2_sol),
             rounded_depth=-1*round((corrected_depth/1))*1) %>%
      rename(lat=Lat,
             long=Long,
             temp=`Temp.oC`,
             oxygen=`Oxygen.uM`,
             fl=Flu,
             pressure=`Pressure..db.`,
             salinity=`Salinity`) %>%
      select(lat,long,temp,
             salinity,sigma_theta,sa,
             ct,corrected_depth,pressure,
             oxygen,fl,o2_sol,o2_sat,rounded_depth) %>%
      group_by(rounded_depth) %>%
      summarize(across(everything(),mean,.names='{.col}')) %>%
      mutate(mld=rounded_depth[min(which(sigma_theta>=(sigma_theta[which(rounded_depth==10)]+0.125)))])
  }
  return(out)
}

## Reading in CTD cast files

ctd_files<-list.files(path='../omz_update_files/omz_stuff',
                      pattern='ETNP1[0-9].*[0-9].csv$',
                      full.names=TRUE)
ctd_summaries<-sapply(ctd_files,function(x) ctd_summarizing_function(read.csv(x)),simplify=FALSE)
ctd_summaries<-list()
for(i in 1:length(ctd_files)){
  ctd_summaries[[i]]<-ctd_summarizing_function(read.csv(ctd_files[i]))
}

## Processing CTD files and doing thermodynamics calculations

ctd_summaries<-lapply(as.list(ctd_files),function(x) ctd_process_2(read.csv(x)))
home_file<-rep(ctd_files,do.call(c,lapply(ctd_summaries,nrow)))
ctd_collapsed<-as.data.frame(cbind(do.call(rbind,ctd_summaries),home_file)) %>%
  mutate(year=gsub('_.*$','',gsub('^.*ETNP','',home_file)),
         station=gsub('.csv','',gsub('^.*_','',home_file)))
station_locations<-sapply(ctd_files,function(x) c(mean(read.csv(x)$Latitude),mean(read.csv(x)$Longitude)))


## Reading in nutrient data
updated_nn<-read.csv('../omz_update_files/etnp_2013_so6.csv') %>%
  mutate(updated_station=paste0('S0',Station),
         year='13')

## Joining nutrient and CTD data
nutrient_and_ctd<-full_join(nutrient_data,ctd_collapsed,by=c('updated_station'='station','year','Depth..m.'='rounded_depth')) %>%
  full_join(updated_nn,by=c('updated_station','year','Depth..m.'='Depth',
                            'Nitrate..uM.','Nitrite..uM.'))

## Selecting example station for figure plotting
environmental_data_plot_df<-nutrient_and_ctd %>%
  filter(updated_station=='S06')

## Nitrogen depth profile

nitrogen_context<-ggplot(drop_na(environmental_data_plot_df,Nitrate..uM.,Nitrite..uM.))+
  geom_point(aes(x=Depth..m.,y=Nitrate..uM.,col='darkblue'),
             size=4)+
  geom_point(aes(x=Depth..m.,y=Nitrite..uM.,col='gold4'),
             size=4)+
  geom_line(aes(x=Depth..m.,y=Nitrate..uM.,col='darkblue'),
            size=1.5)+
  geom_line(aes(x=Depth..m.,y=Nitrite..uM.,col='gold4'),
            size=1.5)+
  scale_color_manual(values=c('darkblue','gold4'),
                     labels=c('Nitrate','Nitrite'),
                     name='')+
  scale_y_continuous(sec.axis=dup_axis(),limits=c(0,24))+
  xlab('Depth [m]')+
  ylab('Nitrogen Species [uM]')+
  coord_flip()+
  scale_x_reverse(lim=c(300,0))+
  theme_bw()+
  theme(text=element_text(size=16),
        legend.position=c(0.8,0.95),
        axis.title.x.bottom=element_blank(),
        axis.text.x.top=element_blank(),
        axis.ticks.x.top=element_blank(),
        legend.background=element_blank(),
        axis.title.y=element_blank())

## Oxygen profile
oxygen_corresponding<-ggplot((filter(nutrient_and_ctd,updated_station=='S06',year==13)))+
  geom_point(aes(x=Depth..m.,y=o2_sat),size=1)+
  coord_flip()+
  scale_x_reverse(lim=c(300,0))+
  scale_y_continuous(sec.axis=dup_axis())+
  theme_bw()+
  xlab('Depth [m]')+
  ylab('Oxygen Saturation %')+
  theme(text=element_text(size=16),
        axis.title.x.bottom=element_blank(),
        axis.text.x.top=element_blank(),
        axis.ticks.x.top=element_blank())

## Reading in metagenomic summary data, including mean per-read GC%
## GET THIS SETTLED
bulk_mat<-read.csv('../gc_OMZ/metagenome_wide_characteristics.csv')
  
gc_profile<-ggplot(bulk_mat)+
  geom_point(aes(x=depth,y=avggc,col=factor(fraction,levels=c('PF','SV'),
                                            labels=c('Particle','Free-Living'))))+
  geom_errorbar(aes(x=depth,ymin=avggc-sqrt(gcvar),ymax=avggc+sqrt(gcvar),
                    col=factor(fraction,levels=c('PF','SV'),
                               labels=c('Particle','Free-Living'))),
                size=3)+
  scale_y_continuous(sec.axis=dup_axis())+
  coord_flip(xlim=c(300,0))+
  ylab('Metagenomic GC%')+
  theme_bw()+
  scale_color_manual(name='Size Fraction',values=c('darkgreen','cyan3'))+
  theme(text=element_text(size=16),
        axis.title.y=element_blank(),
        axis.title.x.bottom=element_blank(),
        axis.text.x.top=element_blank(),
        axis.ticks.x.top=element_blank(),
        legend.position=c(0.775,0.94),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.background=element_blank())

## Conduct model analysis to demonstrate different in bulk GC between size
## fractions
bulk_gc<-lme(fixed=avggc~fraction+log10(nreads),random=~1|depth,
             data=bulk_mat,
             cor=corExp())
## Check out model statistics
summary(bulk_gc)

## Make supplemental table 1
sjPlot::tab_model(bulk_gc,
                  dv.labels='Average GC%',
                  pred.labels=c('Intercept',
                                'Free-Living Size Fraction',
                                'Log10 Sequencing Depth'),
                  file='supplemental_table1.html')

## Compose completed figure 1 and save
full_fig<-map_element/(oxygen_corresponding+nitrogen_context+gc_profile)

ggsave('figure1_full_v2.pdf',full_fig,
       device='pdf',scale=2)


