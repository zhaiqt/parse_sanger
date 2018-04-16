import re
import translator

def read_fasta(inputfilename):
    inputfile = open(inputfilename, 'r')
    dict={}
    name=''
    seq=''
    for line in inputfile:
        line=line.strip('\n')
        if line.startswith('>'):
            if name and seq:
                seq=seq.replace(' ','')
                seq=seq.upper()
                protein=translator.choose_translation(seq)
                dict[name] = {'DNA': seq, 'PRO': protein}
            line=line.lstrip('>')
            #line = re.split(r'\D', line)
            #name = line[0]
            name =line

            seq=''
        elif line:
            seq +=line

    seq=seq.replace(' ','')
    seq=seq.upper()
    print seq
    protein=translator.choose_translation(seq)
    dict[name] = {'DNA': seq, 'PRO': protein}

    #print dict
    print "%d fasta sequences are imported." % len(dict)
    return dict


#########
'''
testname='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/results/LC_consensus_seqs.fasta'
read_fasta(testname)
'''
