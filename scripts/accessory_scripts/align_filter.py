import pysam
import sys
infile=sys.argv[1]
outfile=infile.split('.sort')[0]
outname=outfile+'.edit.bam'
sumname=outfile+'.sum.out'
samfile=pysam.AlignmentFile(infile,"rb")
outfile=pysam.AlignmentFile(outname,"w",template=samfile)
out_summary=open(sumname,"w")
idlist=list()
balist=list()
totalbases=list()
gene_names=list()
for read in samfile.fetch():
	edits=float(read.get_tag("NM"))
	alength=float(read.query_alignment_length)
	qlength=float(read.query_length)
	ref=read.reference_name
	identity=1-edits/alength
	align_percent=alength/qlength
	totalbases.append(qlength)
	if (identity>=0.95) & (alength>60):
		gene_names.append(ref)
		idlist.append(identity)
		balist.append(alength)
		outfile.write(read)
		out_summary.writelines(','.join([str(ref),str(identity),str(alength)])+'\n')
print('Number passed reads: '+str(len(idlist)))
print('Total bases queried: '+str(sum(totalbases)))
print('Percent bases passed and aligned: '+str(sum(balist)/sum(totalbases)))