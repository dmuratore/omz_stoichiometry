## The purpose of this script is to implement our alignment processing pipeline
## for converting metagenome mappings against OMZ gene catalog to coverages, est 
## copy numbers

/usr/local/pacerepov1/samtools/1.0/gcc-4.9.0/bin/samtools view -b -o ${1//.sam/.bam} ${1}
/usr/local/pacerepov1/samtools/1.0/gcc-4.9.0/bin/samtools sort -O BAM -o ${1//.sam/}'.sort.bam' ${1//.sam/.bam} -t 4
/usr/local/pacerepov1/samtools/1.0/gcc-4.9.0/bin/samtools index ${1//.sam/}'.sort.bam' -@ 4
/usr/local/pacerepov1/samtools/1.0/gcc-4.9.0/bin/samtools stats ${1//.sam/}'.sort.bam' > ${1//.sam/.orig.stats.out}
python align_filter.py ${1//.sam/}'.sort.bam'
/usr/local/pacerepov1/samtools/1.0/gcc-4.9.0/bin/samtools sort -O BAM -o ${1//.sam/.sort.edit.bam} ${1//.sam/.edit.bam}
/usr/local/pacerepov1/samtools/1.0/gcc-4.9.0/bin/samtools index ${1//.sam/.edit.bam}
/usr/local/pacerepov1/samtools/1.0/gcc-4.9.0/bin/samtools idxstats ${1//.sam/.edit.bam} > ${1//.sam/_edited_idx.out}
Rscript --vanilla working_with_coverage.R ${1//.sam/_edited_idx.out} ${1//.sam/.sum.out} ${2}