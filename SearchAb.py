import argparse
import os
import translator
import WriteFast
import ReadIgBlastn
import AnnotateProtein
import ReadFasta


parser = argparse.ArgumentParser(prog='process fasta sequences, translate and annotate CDR and germline',
                                 description="python SearchAb.py -s mouse -c h", epilog='')
parser.add_argument('-s', '--species', help='mouse, rabbit or human', default="mouse")
parser.add_argument('-c', '--chain', help="H" or "L", default="H")
parser.add_argument('-k', '--keylist',nargs='+', help="DNA,FV,CDR1,CDR2,CDR3, GERMLINE", default="DNAlen GERMLINE-V CDR3-DNA RID")
parser.add_argument ('-i',"--inputfasta",type=str, help="input files directory", default='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/results/LC_consensus_seqs.fasta',action='store')
parser.add_argument ('-o',"--outputpath",type=str, help="input files directory", default='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/results')

args = parser.parse_args()


########### read fasta file and convert it into dict##########
AbDict = ReadFasta.read_fasta(args.inputfasta)

########### run Igblastn ##########
os.system("python /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/WrapIgBlastn.py -s %s -i %s" %( args.species, args.inputfasta))
igblastnFilename=args.inputfasta.rstrip('.fasta')+".igblastn"

# #############extract results from Igblastn, the results are returned as dictionary {name: }
#'M00680:164:000000000-AN2N4:1:2119                                                     :11126:25130': {'GERMLINE-J': 'VK', 'FR1head_pos': 1, 'CDR3head_pos': 277, 'FR2tail_pos': 159, 'CDR2tail_pos': 168, 'CDR1head_pos': 7                                                     9, 'GERMLINE-D': 'JK1', 'CDR1tail_pos': 108, 'FR2head_pos': 109, 'CDR2head_pos': 160, 'FR1tail_pos': 78, 'GERMLINE-V': '21-4', 'CDR3t                                                     ail_pos': 296},
foo=ReadIgBlastn.ReadIgBlastn(igblastnFilename)
foo.readIgBlastn()
#print igblastn_results
for key in foo._dict:
    #print key
    #print igblastn_results[key]
    AbDict[key].update(foo._dict[key])

##### Anotate Ab Protein 1000 sequences using PWM      #######
foo = AnnotateProtein.AnnotateProtein(AbDict,args.species,args.chain)
foo.AnnotateDict()

'''
print AbDict
for a,b in AbDict.iteritems():
    print a
    print b
'''



######## write all the information of AbDict into all.txt######
### keyList=["GERMLINE-V","GERMLINE-D","GERMLINE-J","PRODUCT","CHAIN","LID","RID","DNA","PRO",'FR1-PRO','CDR1-PRO','FR2-PRO','CDR2-PRO',"FR3-PRO",'CDR3-PRO','FR4-PRO','FR1-DNA','CDR1-DNA','FR2-DNA','CDR2-DNA',"FR3-DNA",'CDR3-DNA','FR4-DNA']
print "write out all the information of AbDict ........"
#WriteFast.writeDict_keys(AbDict,args.outputpath,prefix,keylist)
WriteFast.writeDict_keys(AbDict,args.outputpath,args.species+'_'+args.chain+'_')


'''
1288
{'FR1tail_pos': 200, 'DNA': 'TTATATATGGTGTTATATCAAAGNAGANTGGGGATGTACATCTGAAAGGCAGGTGGAGCAAGATGGAATCACAGACTCAGGTCCTCATCTCCTTGCTGTTCTGGGTATCTGGTACCTGTGGGGACATTGTGATGACACAGTCTCCATCCTCCCTGAGTGTGTCAGCAGGAGAGAAGGTCACTATGAGCTGCAAGTCCAGTCAGAGTCTGTTAAACAGTGGAAATCAAAAGAACTACTTGGCCTGGTACCAGCAGAAACCAGGGCAGCCTCCTAAACTGTTGATCTACGGGGCATCCACTAGGGAATCTGGGGTCCCTGATCGCTTCACAGGCAGTGGATCTGGAACCGATTTCACTCTTACCATCAGCAGTGTGCAGGCTGAAGACCTGGCAGTTTATTACTGTCAGAATGATCATAGTTTTCCGTACACGTTCGGAGGGGGGACCAAGCTGGAAATTAAACGGGCTGATGCTGCACCAACTGTATCCATCTTCCCACCATCCAGTGAGCAGTTAACATCTGGAGGTGCCTCAGTCGTGTGCTTCTGAACAACTCTACCCCNAAA', 'CDR3head_pos': 405, 'CDR2tail_pos': 296, 'FR1-PRO': 'DIVMTQSPSSLSVSAGEKVTMSCKSS', 'CDR2-PRO': 'GAS', 'CDR3-DNA': 'CAGAATGATCATAGTTTTCC', 'CDR1tail_pos': 236, 'FR3-DNA': 'CATCCACTAGGGAATCTGGGGTCCCTGATCGCTTCACAGGCAGTGGATCTGGAACCGATTTCACTCTTACCATCAGCAGTGTGCAGGCTGAAGACCTGGCAGTTTATTACT', 'FR4-PRO': '', 'FR1-DNA': 'GACATTGTGATGACACAGTCTCCATCCTCCCTGAGTGTGTCAGCAGGAGAGAAGGTCACTATGAGCTGCAAGTCCAGT', 'FR2-PRO': 'LAWYQQKPGQPPKLLIY', 'FR3tail_pos': 404, 'CDR1-DNA': 'CAGAGTCTGTTAAACAGTGGAAATCAAAAGAACTAC', 'CDR1-PRO': 'QSLLNSGNQKNY', 'GERMLINE-J': 'VK', 'FR1head_pos': 123, 'CDR3tail_pos': 424, 'FR4-DNA': '', 'FR2tail_pos': 287, 'CDR2-DNA': 'GGGGCATCC', 'PRO': 'IYGVISKxxGDVHLKGRWSKMESQTQVLISLLFWVSGTCGDIVMTQSPSSLSVSAGEKVTMSCKSSQSLLNSGNQKNYLAWYQQKPGQPPKLLIYGASTRESGVPDRFTGSGSGTDFTLTISSVQAEDLAVYYCQNDHSFPYTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCF*TTLPx', 'CDR1head_pos': 201, 'GERMLINE-D': 'JK2', 'FR3head_pos': 297, 'FR2head_pos': 237, 'CDR2head_pos': 288, 'CDR3-PRO': 'QNDHSF', 'FV-DNA': 'GACATTGTGATGACACAGTCTCCATCCTCCCTGAGTGTGTCAGCAGGAGAGAAGGTCACTATGAGCTGCAAGTCCAGTCAGAGTCTGTTAAACAGTGGAAATCAAAAGAACTACTTGGCCTGGTACCAGCAGAAACCAGGGCAGCCTCCTAAACTGTTGATCTACGGGGCATCCCATCCACTAGGGAATCTGGGGTCCCTGATCGCTTCACAGGCAGTGGATCTGGAACCGATTTCACTCTTACCATCAGCAGTGTGCAGGCTGAAGACCTGGCAGTTTATTACTCAGAATGATCATAGTTTTCC', 'FV-PRO': 'DIVMTQSPSSLSVSAGEKVTMSCKSSQSLLNSGNQKNYLAWYQQKPGQPPKLLIYGASSTRESGVPDRFTGSGSGTDFTLTISSVQAEDLAVYYCQNDHSF', 'FR3-PRO': 'STRESGVPDRFTGSGSGTDFTLTISSVQAEDLAVYYC', 'GERMLINE-V': '8-28', 'FR2-DNA': 'TTGGCCTGGTACCAGCAGAAACCAGGGCAGCCTCCTAAACTGTTGATCTAC'}
'''
