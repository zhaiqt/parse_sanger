python Sanger2Fastq -i inputfile.zip -o outputpath

python ParseSanger.py -i /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/raw/TIM1_arm3_HC.zip

python ParseSanger.py -i /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/raw/TIM1_arm3_LC.zip


python SearchAb.py

python SearchAb.py -i /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/results/LC_consensus_seqs.fasta -o /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/results -c L -s mouse

python SearchAb.py -i /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/results/HC_consensus_seqs.fasta  -o /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/parse_sanger/results -c H -s mouse

HC_consensus_seqs.fasta
