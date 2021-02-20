import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time

def combination_signalp_info(signalp_result, reference, direction_path):
    signalp_result = signalp_result + '_summary.signalp5'
    with open(signalp_result, 'r') as file:
        signalp_result_infos = file.readlines()
    file.close()

    record_dic = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))

    n_fasta_seqs = []
    for signalp_result_info in signalp_result_infos[2:]:
        id = signalp_result_info.split()[0]
        tag = signalp_result_info.split()[1]
        if tag != 'OTHER':
            if id in record_dic.keys():
                SP = 'SP(Sec/SPI)=' + signalp_result_info.split()[2]
                TAT = 'TAT(Tat/SPI)=' + signalp_result_info.split()[3]
                LIPO = 'LIPO(Sec/SPII)=' + signalp_result_info.split()[4]
                OTHER = 'OTHER=' + signalp_result_info.split()[5]
                try:
                    CS_Position = signalp_result_info.split('CS pos: ')[1].split('.')[0]
                    if '?' in CS_Position:
                        CS_Position = '000'
                except IndexError:
                    CS_Position = '000'
                CS_Position =  'CS(Position)=' + CS_Position
                position_info = record_dic[id].description.split('[')[1].split(']')[0]
                rec = SeqRecord(
                    Seq(str(record_dic[id].seq)),
                    id=record_dic[id].id,
                    name=record_dic[id].name,
                    description= position_info + '|signalp|' + SP + '|' + TAT + '|' + LIPO + '|' + OTHER + '|' + CS_Position,
                )
                n_fasta_seqs.append(rec)
    SeqIO.write(n_fasta_seqs, direction_path + '/signalp_info_combined.faa', "fasta")


def my_signalp_by_getorf(getorf_result, signalp_org, signalp_format, signalp_result, direction_path):
    signalp_command = 'signalp -fasta ' + getorf_result + ' ' + '-org' + ' ' + signalp_org + ' ' + '-format' + ' ' + signalp_format + ' ' + '-prefix' + ' ' + signalp_result
    os.system(signalp_command)
    combination_signalp_info(signalp_result=signalp_result, reference=getorf_result, direction_path=direction_path)


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

    start_signalp = time.time()
    my_signalp_by_getorf(args.getorf_result, signalp_org='gram-', signalp_format='short',
                         signalp_result=args.direction_path + "/signalp_result", direction_path=args.direction_path)
    done_signalp = time.time()
    elapsed = done_signalp - start_signalp
    print("signalp run: " + str(elapsed / 60.0))

if __name__ == '__main__':
    main()