import zipfile
import os
import TrimEnds
import Merge2FastQ
import ReadSanger
import MergeMultiFastQ
import logging
import csv
import argparse


#python IDAP384screen.py -i /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/IDAP384screen/IDAP_5_star/rawdata/ensingle.csv.zip -c /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/IDAP384screen/IDAP_5_star/rawdata/barcode_Table.csv -o /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/IDAP384screen/IDAP_5_star/result
parser = argparse.ArgumentParser( prog='SangerPairSeqContig',description="Read paired-ends sanger sequneces, convert to fastq, and merge the Forward and the reverse reads, -> 1) fasta contains the assembled sequences, 2) fasta with non-assembled sequences", epilog='python Sanger2Fastq -i inputfile.zip  -o outputpath')
parser.add_argument ('-i','--inputzip',help='Input zipfile', default='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/raw/Archive.zip')
parser.add_argument('-o', '--outputpath',help='outputpath for output',default='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/results')
parser.print_help()
args=parser.parse_args()
print args

log_file = args.outputpath+'/runlog.txt'
log_level = logging.DEBUG
logging.basicConfig(filename=log_file, level=log_level, format='%(asctime)s %(message)s')


####unzip the inputfile to the output filepath#####
abi_path = ReadSanger.extract_zip(args.inputzip,args.outputpath)

#generate directories for the All, F , R fastq files
localinfilename = os.path.basename (abi_path)
localinfilename = localinfilename.strip('.zip')

#local_abi= localinfilename+"_ABI"
#outfilepath_abi = os.path.join(args.outputpath,local_abi)

local_rawfastq = localinfilename+"_raw_fastq"
outfilepath_rawfastq = os.path.join(args.outputpath,local_rawfastq)

local_parsed_fastq = localinfilename+"_parsed_fastq"
outfilepath_parsed_fastq  = os.path.join(args.outputpath,local_parsed_fastq  )

local_group_fasta = localinfilename+"_group_fasta"
outfilepath_group_fasta   = os.path.join(args.outputpath,local_group_fasta  )

local_group_trimed_fasta=localinfilename+"_group_trim_fasta"
outfilepath_group_trimed_fasta = os.path.join(args.outputpath,local_group_trimed_fasta  )

print "1)rawfastq" + " " +outfilepath_rawfastq
print '2)parsed_fastq ' + " " +outfilepath_parsed_fastq

for dirname in [outfilepath_rawfastq,outfilepath_parsed_fastq,outfilepath_group_fasta,outfilepath_group_trimed_fasta  ]:
    try:
        os.makedirs(dirname)
    except OSError:
        if os.path.exists(dirname):
            # We are nearly safe
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise


#####  convert Abi/ab1 to a fastq file  ############
ReadSanger.read_sanger(abi_path,outfilepath_rawfastq)

######## combine all the fastq to a single fastq, And change the name to single dict######
raw_single_fastq_name = os.path.join(args.outputpath,"raw_single.fastq")
all_fastq_dict =ReadSanger.concentrate2single (outfilepath_rawfastq, raw_single_fastq_name)
print 'length of all fastq_dict :' + str(len(all_fastq_dict))


###########
# write untrimed fasta
ReadSanger.write_fasta_bygroup(all_fastq_dict,outfilepath_group_fasta,False)
# write trimed fasta
ReadSanger.write_fasta_bygroup(all_fastq_dict,outfilepath_group_trimed_fasta,True)
# write gouped untrimed fastQ, same direction
ReadSanger.write_fastq_bygroup(all_fastq_dict,outfilepath_parsed_fastq,False)

'''
(All_fastq_path,All_fasta_dict)=MergeMultiFastQ.read_FastQs_dict(outfilepath_All,args.outputpath)

#print len(All_fasta_dict)

All_Fasta_dict = MergeMultiFastQ.Contig_2Fastq_Fasta(All_fasta_dict, args.outputpath)
#print All_Fasta_dict

####### write out the fasta
contig_fasta_name = os.path.join(args.outputpath,'all.fasta')
MergeMultiFastQ.writeFasta (All_Fasta_dict, contig_fasta_name )
'''
