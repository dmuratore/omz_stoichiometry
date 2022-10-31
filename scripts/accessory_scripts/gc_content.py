#The purpose of this script is to calculate the GC content and
#Variation in GC content for each MG sample
#Import necessary packages
import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import csv
print('imports conducted')
#Load input file
in_file=sys.argv[1]
out_dir=sys.argv[2]
in_no_ext=in_file[0:-3]
print(in_file+' loaded')
print(out_dir+' selected')

#Make a variance function because I guess there isn't one
##And also a mean one
def mean(x):
	result=sum(x)/float(len(x))
	return result
def gc_variance(x):
	avg=mean(x)
	result=[(i-avg)**2 for i in x ]
	result=mean(result)
	return result
print('functions defined')

#Initialize vectors for calculating summary statistics
read_gc=[]
read_length=[]
read_gc_percent=[]

#Calculate results
print('Beginning processing')
for seq_record in SeqIO.parse(in_file,"fastq"):
	read_gc_percent.append(GC(seq_record.seq))
	read_length.append(float(len(seq_record)))
	read_gc.append(read_gc_percent[-1]*read_length[-1])
print('Processing complete')

#Do summary statistics
print('Generating Summary Report')
total_gc=sum(read_gc)/sum(read_length)
mean_gc=mean(read_gc_percent)
proportion_gc=[x/100 for x in read_gc_percent]
var_gc=gc_variance(proportion_gc)
m_read_length=mean(read_length)
n_reads=len(read_gc)
summary_report=[in_no_ext,total_gc,mean_gc,var_gc,m_read_length,n_reads]
filename=out_dir+in_no_ext+'_summary.csv'
with open(filename,'w') as f:
	settings=csv.writer(f,delimiter=',')
	settings.writerow(summary_report)
print('Report Complete')



