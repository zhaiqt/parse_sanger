
from Bio import SeqIO
import zipfile
import os
import re
import translator
import TrimEnds

##################
def extract_zip(infile,outfilePath):
    if infile.endswith(".zip"):
        with zipfile.ZipFile(infile, 'r') as z:
            print "zipfile zipfile"
            localfilename = os.path.basename(z.filename)
            abi_path = os.path.join(outfilePath,localfilename.rsplit('.zip')[0])
            z.extractall(abi_path)
        print abi_path
    return abi_path



###############
def read_sanger(infilepath,outfilepath):
    count_infile = 0
    print infilepath
    for filename in os.listdir(infilepath):
        #print "infilepath :" +infilepath
        #print filename
        if filename.endswith('.ab1'):
            abi_filename=os.path.join(infilepath,filename)
            fastq_filename= os.path.join(outfilepath,filename.replace('.ab1','.fastq'))
            #print abi_filename
            ##print fastq_filename
            SeqIO.convert(abi_filename,'abi',fastq_filename,'fastq')
            '''
            if filename.endswith("QB5505.ab1"):
                F_fastq_filename =os.path.join(outfilepath_F,filename.replace('.ab1','.fastq'))
                SeqIO.convert(abi_filename,'abi',F_fastq_filename,'fastq')
            else:
                R_fastq_filename =os.path.join(outfilepath_R,filename.replace('.ab1','.fastq'))
                SeqIO.convert(abi_filename,'abi',R_fastq_filename,'fastq')
            '''
            count_infile += 1
    print "There are total %d sequences" %count_infile

    return

##########



############################################
def concentrate2single(infilepath,outputfilename):  #also reverse translate
    #combinedfile_name = os.path.basename(infilepath)
    raw_single_fasta_name =outputfilename.rstrip('.fastq')+'.fasta'
    raw_single_fasta_file = open(raw_single_fasta_name, 'wb')

    all_fastq = {}
    count_file =0
    #combinedfile_name =os.path.join(outputfilename,combinedfile_name+'.fastq')
    combinedfile_name = outputfilename
    print "concentrate2single name is" +combinedfile_name
    combinedfile = open(combinedfile_name , 'wb')

    print outputfilename
    for filename in os.listdir(infilepath):
        single_fastq = []
        flag_reverse = False

        if filename.endswith('.fastq'):
            seq_name = filename.split(';')[0]
            with open(os.path.join(infilepath,filename)) as f:
                for row in f:
                    row = row.strip('\n')
                    if row.startswith("@"):
                        fasta_output = '\n>'+row.lstrip('@')+'\n'
                        # row = row.split(";")
                        if "QB6179" in row or "QB5506" in row or 'QB6178' in row or 'Rev' in row or 'rev' in row:
                            flag_reverse = True
                        IDs=re.search('\D(\d{4})\D',row)  ###### Extract 4 numbers , return
                        ID=IDs.group(0)  #####
                        #row= row[0]+"\n"
                        #print ID
                        row = "@"+ID
                        #ID=str(ID) +'\n'
                    combinedfile.write(row+'\n')
                    single_fastq.append(row.lstrip('@'))

                single_fastq[-1] = covert_Qscore(single_fastq[-1])
                fasta_output +=single_fastq[1] + '\n'
                raw_single_fasta_file.write(fasta_output)


                if flag_reverse == True:
                    try:
                        single_fastq =translator.reverse_complement_FastQ(single_fastq)
                    except:
                        print single_fastq
                all_fastq[ID] = all_fastq.get(ID, [])
                all_fastq[ID].append(single_fastq)

                f.close()
    combinedfile.close()
    print "All the fastqs were concentrated into a single file ------%s. And the ID was extracted." +combinedfile_name
    return all_fastq



########
def covert_Qscore(String):
	Q_score=[]
	for i in String:
		Q_score.append(int((ord(i))-33))

        #Q_score ([ord(x)-33 for x in String])
	return Q_score

########
def reverse_convert_Qscore(list):
    Q_string=''
    for i in list:
        Q_string += chr(i+33)
    return Q_string



##############
def write_fasta_bygroup(indict,outputpath,flag_trim):

    for ID,fastqs in indict.iteritems():
        fasta_outname = os.path.join(outputpath,ID)
        fasta_outfile = open(fasta_outname,'w')

        for fastq in fastqs:
            if flag_trim == True:

                object1 = TrimEnds.TrimEnds(fastq)
                object1.trim5End(minimum_quality=20)
                fastq = object1.output_trimed_fastq()

            output='>'+fastq[0]+'\n'+fastq[1]+'\n'
            fasta_outfile.write(output)
        fasta_outfile.close()
    return
####################

##############
def write_fastq_bygroup(indict,outputpath,flag_trim):

    for ID,fastqs in indict.iteritems():
        fastq_outname = os.path.join(outputpath,ID)
        fastq_outfile = open(fastq_outname,'w')
        i=0
        for fastq in fastqs:
            if flag_trim == True:

                object1 = TrimEnds.TrimEnds(fastq)
                object1.trim5End(minimum_quality=20)
                fastq = object1.output_trimed_fastq()

            output='@'+fastq[0]+'_'+str(i)+'\n'+fastq[1]+'\n'+fastq[2]+'\n'+reverse_convert_Qscore(fastq[3])+'\n'
            fastq_outfile.write(output)
            i +=1

        fastq_outfile.close()
    return
