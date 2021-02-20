import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time

def combination_tmhmm_info(tmhmm_result, reference, direction_path):
    with open(tmhmm_result, 'r') as file:
       tmhmm_result_infos = file.readlines()
    file.close()
    record_dic = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))

    n_fasta_seqs = []
    for tmhmm_result_info in tmhmm_result_infos:
        id = tmhmm_result_info.split()[0]
        PredHel_value = tmhmm_result_info.split()[4].split("PredHel=")[1]
        if int(PredHel_value) > 0:
            if id in record_dic.keys():
                length = tmhmm_result_info.split()[1]
                ExpAA = tmhmm_result_info.split()[2]
                First60 = tmhmm_result_info.split()[3]
                PredHel = tmhmm_result_info.split()[4]
                position_info = record_dic[id].description.split('[')[1].split(']')[0]
                rec = SeqRecord(
                    Seq(str(record_dic[id].seq)),
                    id=record_dic[id].id,
                    name=record_dic[id].name,
                    description=position_info+'|tmhmm|'+length+'|'+ExpAA+'|'+First60+'|'+PredHel,
                )
                n_fasta_seqs.append(rec)
    SeqIO.write(n_fasta_seqs, direction_path + '/tmhmm_info_combined.faa', 'fasta')


def my_tmhmm_getorf(getorf_result, direction_path):
    tmhmm_result = direction_path + '/tmhmm_result.txt'
    tmhmm_command = './tmhmm-2.0c/bin/tmhmm ' + getorf_result + ' -short > ' + tmhmm_result
    print("tmhmm coming.......")
    try:
        os.system(tmhmm_command)
        os.system('rm -rf TMHMM_*')
    except:
        print("tmhmm error")
    combination_tmhmm_info(tmhmm_result=tmhmm_result, reference=getorf_result, direction_path=direction_path)


def my_getorf(getorf_in, getorf_result, getorf_table, getorf_minsize, getorf_maxsize):
    getorf_command = 'getorf -sequence ' + getorf_in + ' ' +  '-outseq' + ' ' + getorf_result + ' ' + '-table' + ' ' + getorf_table
    getorf_command = getorf_command + ' ' + '-minsize' + ' ' + getorf_minsize + ' ' + '-maxsize' + ' ' + getorf_maxsize
    os.system(getorf_command)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--direction_path", type=str, required=True)
    parser.add_argument("-gi", "--getorf_in", type=str, required=False, default= 'GCF_003018455.1_ASM301845v1_genomic.fna')
    parser.add_argument("-gr", "--getorf_result", type=str, required=False, default= 'GCF_003018455.1_ASM301845v1_genomic.ORF.15-50aa.faa')
    parser.add_argument("-gt", "--getorf_table", type=str, required=False, default='1')
    parser.add_argument("-gmi", "--getorf_minsize", type=str, required=False, default='45')
    parser.add_argument("-gma", "--getorf_maxsize", type=str, required=False, default='150')
    args = parser.parse_args()

    path = args.direction_path
    try:
        os.mkdir(path)
    except:
        pass

    start_getorf = time.time()
    my_getorf(args.getorf_in, args.getorf_result, args.getorf_table, args.getorf_minsize, args.getorf_maxsize)
    done_getorf = time.time()
    elapsed = done_getorf - start_getorf
    print("getorf run: " + str(elapsed/60.0))


    args.getorf_result = path + '/' + args.getorf_result
    start_tmhmm = time.time()
    my_tmhmm_getorf(args.getorf_result, args.direction_path)
    done_tmhmm = time.time()
    elapsed = done_tmhmm - start_tmhmm
    print("tmhmm run: " + str(elapsed / 60.0))

if __name__ == '__main__':
    main()