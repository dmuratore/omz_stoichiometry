## The purpose of this script is to find the MW/N-/C-ARSC for the predicted genes
## from the 13 14 OMZ data
import sys
sys.path.append('/nv/hp10/dmuratore3/data/packages/gene-characteristics-master/')
import codon_table as ct
from Bio import SeqIO, SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis
in_seqs=sys.argv[1]
out_stats=sys.argv[2]
open_inseqs=open(in_seqs,"rU")
gene_set=SeqIO.parse(open_inseqs,"fasta")
errorfile_handle=open('/nv/hp10/dmuratore3/data/job_outputs/stoich_errors.out','w')
outfile=open(out_stats,'w')
aa_list=['A','C','E','D','G','F','I','H','K','M','L','N','Q','P','S','R','T','W','V','Y']
outlist=['gene.id','gene_length','gene_GC','Codon_usage','N_ARSC','C_ARSC','N_C_rat','MW','choice_rank']+aa_list
outfile.writelines(','.join(outlist)+'\n')

timer=0
while True:
	try:
		gene=gene_set.next()
		timer+=1

		if len(gene.seq)<100:
			error=gene.id.strip()+' is <100 bp'
			errorfile_handle.writelines(error)
			continue

		if ct.scan_for_stop_codons(gene.seq)=='True':
			error=gene.id.strip()+' has internal stop codons'
			errorfile_handle.writelines(error)
			continue

		gene_GC=str(SeqUtils.GC(gene.seq[0:-3]))
		N_ARSC,MW,C_ARSC,N_C_rat=ct.ARSC_MW_from_nucleotides(ct.screen_out_ambiguous_codons(gene.seq))
		Codon_usage,choice_rank=ct.calculate_SCU(gene,errorfile_handle)
		protein_seq=str(gene.seq.translate())
		protein_seq=ProteinAnalysis(protein_seq)
		aa_content=protein_seq.count_amino_acids()
		aa_content=[str(x) for x in aa_content.values()]
		if len(gene.seq)<300:
			Codon_usage='Too short'
		gene_length=str(len(gene.seq[0:-3]))
		profile=[gene.id,gene_length,gene_GC,Codon_usage,N_ARSC,C_ARSC,N_C_rat,MW,choice_rank]+aa_content
		outfile.writelines(','.join(profile)+'\n')

	except StopIteration as e:
		break
