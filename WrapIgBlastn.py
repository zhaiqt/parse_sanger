# Purpse: To execute IgBlastn without having to type the whole command each time
#
# Functionality: Replace certain variable parts of the commands with user input and execute IgBlast
# Note: This script has to be withing the IgBlast folder so the folder structure has to follow the IgBlast structure
#
# Usage: python WrapIgBlastn.py -s species -i Input_filename

#from argparse import ArgumentParser
from argparse import ArgumentParser
import os

os.chdir('/Users/zhaiqi1/Documents/Novartis/my_code/NGS')
print "new diectory:"+ os.getcwd()

# Parsing arguments
parser=ArgumentParser()
parser.add_argument("-s","--species",default='mouse',help='Enter the species you are interested. Example: human, mouse, rabbit. Default: mouse')
parser.add_argument("-i","--input",default="/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/results/LC_consensus_seqs.fasta", help='Input FASTA file')
args=parser.parse_args()
#from IPython import embed
#embed()




# Executing IgBlast command
#os.system("./igblastp -germline_db_V database/%s_gl_V  -domain_system imgt -num_alignments_V 1 -outfmt 3 -num_threads 2 -query %s -organism %s > %s-igblast.txt" %(species,species,input,input))
command = "igblastn -germline_db_V  /Users/zhaiqi1/Documents/Novartis/my_code/NGS/database/%s_gl_V -germline_db_D /Users/zhaiqi1/Documents/Novartis/my_code/NGS/database/%s_gl_D -germline_db_J /Users/zhaiqi1/Documents/Novartis/my_code/NGS/database/%s_gl_J -domain_system imgt -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -outfmt 7 -auxiliary_data /Users/zhaiqi1/Documents/Novartis/my_code/NGS/optional_file/%s_gl.aux -query %s -organism %s > %s.igblastn" %(args.species,args.species,args.species,args.species,args.input,args.species,args.input.rstrip('.fasta'))
os.system(command)
print command
#os.system("/Users/zhaiqi1/Documents/Novartis/my_code/NGS/1.4.0/ncbi-igblast-1.4.0/bin/igblastn -germline_db_V  /home/zhaiqi1/NGS/IgBlast/database/%s_gl_V -germline_db_D /home/zhaiqi1/NGS/IgBlast/database/%s_gl_D -germline_db_J /home/zhaiqi1/NGS/IgBlast/database/%s_gl_J -domain_system imgt -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -outfmt 7 -auxiliary_data /home/zhaiqi1/NGS/IgBlast/optional_file/%s_gl.aux -query %s -organism %s > %s.igblastn" %(args.species,args.species,args.species,args.species,args.input,args.species,args.input.rstrip('.fasta')))
print "##################"
print "IgBlastn has been performed. the results are in %s.igblastn" % args.input.rstrip('.fasta')
print "################\n"
