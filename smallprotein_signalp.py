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
    cutted_fasta_seqs = []
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
                    origi_CS_Position = signalp_result_info.split('CS pos: ')[1].split('.')[0]
                    if '?' in origi_CS_Position:
                        check_CS_Position = '000'
                    else:
                        check_CS_Position = origi_CS_Position
                except IndexError:
                    check_CS_Position = '000'
                CS_Position =  'CS(Position)=' + check_CS_Position
                position_info = record_dic[id].description.split('[')[1].split(']')[0]
                rec = SeqRecord(
                    Seq(str(record_dic[id].seq)),
                    id=record_dic[id].id,
                    name=record_dic[id].name,
                    description= position_info + '|signalp|' + SP + '|' + TAT + '|' + LIPO + '|' + OTHER + '|' + CS_Position,
                )
                n_fasta_seqs.append(rec)

                if check_CS_Position != '000':
                    pos = int(check_CS_Position.split('-')[0])
                    original_seq = str(record_dic[id].seq)
                    cutted_seq = original_seq[pos:]
                    cutted_rec = SeqRecord(
                        Seq(cutted_seq),
                        id=record_dic[id].id,
                        name=record_dic[id].name,
                        description=position_info + '|signalp|' + SP + '|' + TAT + '|' + LIPO + '|' + OTHER + '|' + CS_Position,
                    )
                    cutted_fasta_seqs.append(cutted_rec)

    SeqIO.write(n_fasta_seqs, direction_path + '/signalp_info_combined.faa', "fasta")
    SeqIO.write(cutted_fasta_seqs, direction_path + '/cut_signalp_info_combined.faa', "fasta")

def my_signalp_by_getorf(getorf_result, signalp_org, signalp_format, signalp_result, direction_path):
    signalp_command = 'signalp -fasta ' + getorf_result + ' ' + '-org' + ' ' + signalp_org + ' ' + '-format' + ' ' + signalp_format + ' ' + '-prefix' + ' ' + signalp_result
    os.system(signalp_command)
    combination_signalp_info(signalp_result=signalp_result, reference=getorf_result, direction_path=direction_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--direction_path", type=str, required=True)
    parser.add_argument("-gr", "--getorf_result", type=str, required=False, default= 'GCF_003018455.1_ASM301845v1_genomic.ORF.15-50aa.faa')
    args = parser.parse_args()

    path = args.direction_path
    try:
        os.mkdir(path)
    except:
        pass

    start_signalp = time.time()
    getorf_result = args.direction_path + '/' + args.getorf_result
    my_signalp_by_getorf(getorf_result, signalp_org='gram-', signalp_format='short',
                         signalp_result=args.direction_path + "/signalp_result", direction_path=args.direction_path)
    done_signalp = time.time()
    elapsed = done_signalp - start_signalp
    print("signalp run: " + str(elapsed / 60.0))

if __name__ == '__main__':
    main()
    # signalp_result = 'result_Feb18/signalp_result'
    # getorf_result = 'result_Feb18/GCF_003018455.1_ASM301845v1_genomic.ORF.15-50aa.faa'
    # direction_path = 'mine'
    # combination_signalp_info(signalp_result=signalp_result, reference=getorf_result, direction_path=direction_path)