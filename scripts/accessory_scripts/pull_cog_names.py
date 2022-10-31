## Let's see if we can get the names of the COGs and turn it into a comprehensible thing
import sys
input_file=sys.argv[1]
input_header=input_file.split(".")[0]
output_name=input_header+".out"
with open(output_name,"w") as out:
	with open(input_file) as f:
		i=1
		for line in f:
			line=line[0:-1]
			outline=line+'\t'+input_header+'_'+str(i)+'\n'
			out.writelines(outline)
			i=i+1