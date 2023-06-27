import sys
import subprocess
import os

os.chdir('/rds/general/project/hda-22-23/live/Summer_projects/dds122/')

# generate kmer counts
k = 21
fastaFile = 'pone.0106545.s006.fasta'
kmerCmd = 'kmc -k%d -fm %s out2 tmp' % (k, fastaFile)
kmcToolsCmd = 'kmc_tools transform out2 dump kmer_count_output21.txt'

try:
    subprocess.check_output(kmerCmd, shell=True)
    subprocess.check_output(kmcToolsCmd, shell=True)
    sys.stdout.write("%s\n" % (str(result)))
except subprocess.CalledProcessError as error:
    sys.stdout.write("%s\n" % (str(error)))

# create dictionaries
kmer_count_output = open('/rds/general/project/hda-22-23/live/Summer_projects/dds122/kmer_count_output21.txt', 'r')
dict_output = {}
count = 0
for line in data.split('\n'):
    if len(line.strip().split('\t')) == 2:
        (sequence_id, count_dict) = line.strip().split('\t')
        sequence_id = sequence_id[1:]
        kmers = dict((k,int(v)) for (k,v) in [d.split(':') for d in count_dict.split(' ')])
        dict_output[sequence_id] = kmers
    else: # if count is on a new line
        if line.startswith('>'):
            start_of_line = line
            count += 1
        else:
            end_of_line = line
            count += 1
        
        if count%2 == 0:
            new_line = start_of_line + end_of_line
            (sequence_id, count_dict) = new_line.strip().split('\t')
            sequence_id = sequence_id[1:]
            kmers = dict((k,int(v)) for (k,v) in [d.split(':') for d in count_dict.split(' ')])
            dict_output[sequence_id] = kmers
