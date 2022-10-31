library(optparse)
library(data.table)
library(svMisc)
option_list<-list(
  make_option('--idx'),
  make_option('--sum'),
  make_option('--out')
)
args=parse_args(OptionParser(option_list = option_list))
filename=gsub('.*/','',args$idx)
filename=gsub('_.*','',filename)
meow<-fread(args$idx,header=FALSE,sep='\t')
print('Read in first file')
read_lengths<-fread(args$sum,header=FALSE,stringsAsFactors = FALSE,sep=',')
print('Read in second file')
genes<-unique(read_lengths$V1)
bases_mapped<-c()
coverages<-c()
read_ids<-read_lengths$V1
read_bases<-read_lengths$V3
meow_ids<-meow$V1
meow_covs<-meow$V2
for(i in 1:length(genes)){
  progress(i,max.value=length(genes))
  gene_lines<-which(read_ids==genes[i])
  bases_mapped[i]<-sum(read_bases[gene_lines])
  coverages[i]<-bases_mapped[i]/meow_covs[which(meow_ids==genes[i])]
}
print('Completed Coverage Calculation')
cogs<-c("COG0012","COG0016","COG0018","COG0172","COG0215","COG0495",
        "COG0525","COG0533","COG0541","COG0552")
cog_covs<-c()
for(i in 1:length(cogs)){
  cog_genes<-grep(cogs[i],genes)
  cog_covs[i]<-sum(coverages[cog_genes])
}
mean_cog_cov<-mean(cog_covs)
copy_number<-coverages/mean_cog_cov
print('Completed gene copy number estimation')
#est_genome_size<-sum(copy_number)/1e3
#write.csv(data.frame(genes,coverages,copy_number),file=args[3])

gene_cat<-fread(args$out,header=TRUE,sep=',')
#gene_cat<-read.csv('~/Desktop/omz_stuff/omz_gene_stats2.csv',header=TRUE)
new_sample_cov<-rep(0,nrow(gene_cat))
new_sample_cn<-rep(0,nrow(gene_cat))
for(i in 1:length(genes)){
  row<-which(gene_cat$gene.id==genes[i])
  new_sample_cov[row]<-coverages[i]
  new_sample_cn[row]<-copy_number[i]
}
print('Completed information reorganization')
updated_gene_cat<-cbind(new_sample_cov,new_sample_cn)
columns<-ncol(updated_gene_cat)
colnames(updated_gene_cat)[columns-1]<-paste0(filename,'_cov')
colnames(updated_gene_cat)[columns]<-paste0(filename,'_cn')
write.csv(updated_gene_cat,file=paste0(filename,'_col.csv'))
print('Calculation complete')
#mapped_reads<-sum(meow[,3])
#unmapped_reads<-sum(meow[,4])
#percent_reads_mapped<-mapped_reads/(unmapped_reads+mapped_reads)
#genes_with_maps<-meow[which(meow[,3]!=0),]
#percent_genes_mapped<-nrow(genes_with_maps)/nrow(meow)
#reads_mapped_per_base<-genes_with_maps[,3]/genes_with_maps[,2]
#meow2<-order(reads_mapped_per_base,decreasing=TRUE)
#meow3<-genes_with_maps[meow2,]
#meow4<-meow3[,3]/meow3[,2]
#cogs<-c("COG0012","COG0016","COG0018","COG0172","COG0215","COG0495",
#        "COG0525","COG0533","COG0541","COG0552")
#maps<-c()
#thits<-c()
#for(i in 1:length(cogs)){
#  hits<-grep(cogs[i],meow3[,1])
#  total_hits<-sum(meow3[hits,3])
#  thits[i]<-total_hits
#  normed_hits<-sum(meow3[hits,3]/meow3[hits,2])
#  maps[i]<-normed_hits
#  print(total_hits)
#  print(normed_hits)
#}
#meow5<-meow4/mean(maps)
#est_g_size<-sum(meow5)/1e3
#rpkm<-meow3[,3]*(1000/meow3[,2])*(1e6/sum(meow3[,3]))
