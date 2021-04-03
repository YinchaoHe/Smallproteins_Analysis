import sys
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
        PredHel_value = tmhmm_result_info.split()[4].split("PredHel=")[1]
        if int(PredHel_value) > 0:
            id = tmhmm_result_info.split()[0]
            pos_hel = tmhmm_result_info.split()[-1]
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
                    description=position_info+'|tmhmm|'+length+'|'+ExpAA+'|'+First60+'|'+PredHel+ '|' +pos_hel,
                )
                n_fasta_seqs.append(rec)
    SeqIO.write(n_fasta_seqs, direction_path + '/tmhmm_info_combined.faa', 'fasta')


def my_tmhmm_getorf(getorf_result, direction_path):
    tmhmm_result = direction_path + '/tmhmm_result.txt'
    tmhmm_command = 'tmhmm ' + getorf_result + ' -short > ' + tmhmm_result
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


def getorf_tmhmm(path, getorf_in, file):
    getorf_result = path + '/' + file + '.ORF.15-50aa.faa'
    getorf_table = '1'
    getorf_minsize = '45'
    getorf_maxsize = '150'

    start_getorf = time.time()
    my_getorf(getorf_in, getorf_result, getorf_table, getorf_minsize, getorf_maxsize)
    done_getorf = time.time()
    elapsed = done_getorf - start_getorf


    #args.getorf_result = path + '/' + args.getorf_result

    my_tmhmm_getorf(getorf_result, path)
    done_tmhmm = time.time()
    elapsed = done_tmhmm - start_getorf
    print("tmhmm run: " + str(elapsed / 60.0))

    command = 'echo ' + getorf_in
    os.system(command)
    command = 'echo' + "totally run: " + str(elapsed / 60.0)
    os.system(command)

def main():
    sec_level_folder_path = sys.argv[1] + '/'
    third_level_folders = os.listdir(sec_level_folder_path)
    for third_level_folder in third_level_folders:
        third_level_folder_path = sec_level_folder_path + third_level_folder + '/'
        forth_level_folders = os.listdir(third_level_folder_path)
        for forth_level_folder in forth_level_folders:
            if 'tmhmm_result.txt' == forth_level_folder:
                command = 'rm ' + third_level_folder_path + forth_level_folder
                os.system(command)
                continue
            forth_level_folder_path = third_level_folder_path + forth_level_folder + '/'
            files = os.listdir(forth_level_folder_path)
            if len(files) < 4:
                file = forth_level_folder + '.faa'
                file_path = forth_level_folder_path + file

                try:
                    getorf_tmhmm(path=forth_level_folder_path, getorf_in=file_path, file=file)
                except:
                    with open('getorf_signalp_error_report.txt', 'a+') as f:
                        f.write(file_path + '\n')
                    f.close()
            else:
                command = 'echo ' + '******DONE******' + forth_level_folder_path
                os.system(command)




if __name__ == '__main__':
    main()
