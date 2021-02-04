import argparse
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_ids(in_file_name, out_id_file_name):
    try:
        os.mkdir('results')
    except:
        pass
    with open(in_file_name , 'r') as file:
        signalp_results = file.readlines()
    file.close()

    target_ids = []
    output_file = open('results/' + out_id_file_name, 'w')
    amount = 0
    for signalp_result in signalp_results[1:]:
        id = signalp_result.split()[0]
        target_ids.append(id)
        output_file.write(id + '\n')
    output_file.close()

def main(in_file_name, out_id_file_name, fasta_file):
    extract_ids(in_file_name, out_id_file_name)
    with open('results/' + out_id_file_name, 'r') as f:
        ids = f.readlines()
    f.close()

    search_ids = []
    for id in ids:
        id = id.split('\n')[0]
        search_ids.append(id)

    records = list(SeqIO.parse("cdhit.result", "fasta"))
    n_fasta_seqs = []
    for record in records:
        if record.id in search_ids:
            print(record)
            rec = SeqRecord(
                Seq(str(record.seq)),
                id = record.id,
                name = record.name,
                description = record.description,
            )
            n_fasta_seqs.append(rec)
    SeqIO.write(n_fasta_seqs, "results/" + fasta_file, "fasta")

def count():
    records = list(SeqIO.parse("cdhit.result", "fasta"))
    print(len(records))

if __name__ == '__main__':
    count()
    exit()
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, required=True)
    parser.add_argument("-id", "--id", type=str, required=True)
    parser.add_argument("-f", "--fasta", type=str, required=True)
    args = parser.parse_args()
    main(in_file_name=args.infile, out_id_file_name=args.id, fasta_file=args.fasta)
