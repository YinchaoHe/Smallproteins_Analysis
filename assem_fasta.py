import argparse
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_ids(path, in_file_name, out_id_file_name):
    try:
        os.mkdir(path + '/txt_processing')
    except:
        pass
    with open(in_file_name , 'r') as file:
        input_results = file.readlines()
    file.close()

    target_ids = []
    output_file = open(path + '/txt_processing/' + out_id_file_name, 'w')
    if 'signalp5' in in_file_name:
        for input_result in input_results[2:]:
            id = input_result.split()[0]
            target_ids.append(id)
            output_file.write(id + '\n')
        output_file.close()
    else:
        for input_result in input_results:
            id = input_result.split()[0]
            target_ids.append(id)
            output_file.write(id + '\n')
        output_file.close()

def seq_len_signalp(signalp_result, cdhit_result):
    signalp_result = signalp_result + '_summary.signalp5'
    with open(signalp_result, 'r') as file:
       input_results = file.readlines()
    file.close()

    records = list(SeqIO.parse(cdhit_result, "fasta"))

    add_length_result = [input_results[0]]
    title = input_results[1].split("\n")[0] + " " * 4 + 'sequence_length' + '\n'
    add_length_result.append(title)
    for input_result in input_results[2:]:
        id = input_result.split()[0]
        for record in records:
            if id == record.id:
                leng_seq = input_result.split("\n")[0] + " " * 4 + str(len(record.seq)) + '\n'
                add_length_result.append(leng_seq)
    with open('add_Seqlength_'+signalp_result, 'w') as f:
        f.writelines(add_length_result)
    f.close()

def filter_signalp(signalp_result):
    signalp_result = signalp_result + '_summary.signalp5'
    with open(signalp_result, 'r') as file:
        input_results = file.readlines()
    file.close()

    signalp_result_name_blocks = signalp_result.split('/')
    path = signalp_result_name_blocks[0]
    filter_signalp_result = 'filtered_' + signalp_result_name_blocks[1]
    output_file = open( path + '/' + filter_signalp_result, 'w')
    output_file.write(input_results[0])
    output_file.write(input_results[1])

    for input_result in input_results[2:]:
        prediction = input_result.split()[1]
        if prediction != 'OTHER':
            output_file.write(input_result)
    output_file.close()


def assemble_fasta(path, signalp_file_name, out_id_file_name, fasta_file, cdhit_result):
    extract_ids(path, signalp_file_name, out_id_file_name)
    with open(path + '/txt_processing/' + out_id_file_name, 'r') as f:
        ids = f.readlines()
    f.close()

    search_ids = []
    for id in ids:
        id = id.split('\n')[0]
        search_ids.append(id)

    records = list(SeqIO.parse(cdhit_result, "fasta"))
    n_fasta_seqs = []
    for record in records:
        if record.id in search_ids:
#            print(record)
            rec = SeqRecord(
                Seq(str(record.seq)),
                id = record.id,
                name = record.name,
                description = record.description,
            )
            n_fasta_seqs.append(rec)
    SeqIO.write(n_fasta_seqs, path + '/' + fasta_file, "fasta")

def count_fasta(path, fasta_file):
    records = list(SeqIO.parse(path + '/' + fasta_file, "fasta"))
    print(len(records))

if __name__ == '__main__':
    seq_len_signalp()
    exit()
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, required=True)
    parser.add_argument("-id", "--id", type=str, required=True)
    parser.add_argument("-f", "--fasta", type=str, required=True)
    args = parser.parse_args()
