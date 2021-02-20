import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
import threading


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

def count_fasta(path, fasta_file):
    records = list(SeqIO.parse(path + '/' + fasta_file, "fasta"))
    print(len(records))

def my_getorf(getorf_in, getorf_result, getorf_table, getorf_minsize, getorf_maxsize):
    getorf_command = 'getorf -sequence ' + getorf_in + ' ' +  '-outseq' + ' ' + getorf_result + ' ' + '-table' + ' ' + getorf_table
    getorf_command = getorf_command + ' ' + '-minsize' + ' ' + getorf_minsize + ' ' + '-maxsize' + ' ' + getorf_maxsize
    os.system(getorf_command)

def my_cdhit(cdhit_input, cdhit_result, cdhit_n, cdhit_p, cdhit_c, cdhit_d, cdhit_M, cdhit_l, cdhit_s, cdhit_aL, cdhit_g):
    cdhit_command = 'cd-hit -i ' + cdhit_input + ' ' + '-o' + ' ' + cdhit_result + ' ' + '-n' + ' ' + cdhit_n + ' ' + '-p' + ' ' + cdhit_p
    cdhit_command = cdhit_command +  ' ' + '-c' + ' ' + cdhit_c + ' ' + '-d' + ' ' + cdhit_d + ' ' + '-M' + ' ' + cdhit_M + ' ' + '-l' + ' ' + cdhit_l + ' ' + '-s' + ' ' + cdhit_s + ' ' + '-aL' + ' ' + cdhit_aL + ' ' + '-g' + ' ' + cdhit_g
    os.system(cdhit_command)


def my_signalp_by_getorf(getorf_result, signalp_org, signalp_format, signalp_result, direction_path):
    signalp_command = 'signalp -fasta ' + getorf_result + ' ' + '-org' + ' ' + signalp_org + ' ' + '-format' + ' ' + signalp_format + ' ' + '-prefix' + ' ' + signalp_result
    os.system(signalp_command)
    combination_signalp_info(signalp_result=signalp_result, reference=getorf_result, direction_path=direction_path)

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


def all_tasks():
    path = 'intermediate'
    try:
        os.mkdir(path)
    except:
        pass
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--direction_path", type=str, required=True)
    parser.add_argument("-g", "--getorf", type=bool, required=False, default=False)
    parser.add_argument("-gi", "--getorf_in", type=str, required=False, default= 'GCF_003018455.1_ASM301845v1_genomic.fna')
    parser.add_argument("-gr", "--getorf_result", type=str, required=False, default= 'GCF_003018455.1_ASM301845v1_genomic.ORF.15-50aa.faa')
    parser.add_argument("-gt", "--getorf_table", type=str, required=False, default='1')
    parser.add_argument("-gmi", "--getorf_minsize", type=str, required=False, default='15')
    parser.add_argument("-gma", "--getorf_maxsize", type=str, required=False, default='50')

    parser.add_argument("-c", "--cdhit", type=bool, required=False, default=False)
    parser.add_argument("-ci", "--cdhit_input", type=str, required=False, default='None')
    parser.add_argument("-cr", "--cdhit_result", type=str, required=False, default='cdhit.result')
    parser.add_argument("-cn", "--cdhit_n", type=str, required=False, default='2')
    parser.add_argument("-cp", "--cdhit_p", type=str, required=False, default='1')
    parser.add_argument("-cc", "--cdhit_c", type=str, required=False, default="0.5")
    parser.add_argument("-cd", "--cdhit_d", type=str, required=False, default="200")
    parser.add_argument("-cM", "--cdhit_M", type=str, required=False, default="50000")
    parser.add_argument("-cl", "--cdhit_l", type=str, required=False, default="5")
    parser.add_argument("-cs", "--cdhit_s", type=str, required=False, default="0.95")
    parser.add_argument("-caL", "--cdhit_aL", type=str, required=False, default="0.95")
    parser.add_argument("-cg", "--cdhit_g", type=str, required=False, default="1")


    parser.add_argument("-s", "--signalp", type=bool, required=False, default=False)
    parser.add_argument("-sg", "--signalp_getorf", type=bool, required=False, default=False)
    parser.add_argument("-si", "--signalp_input", type=str, required=False, default='None')
    parser.add_argument("-so", "--signalp_org", type=str, required=False, default="gram-" )
    parser.add_argument("-sf", "--signalp_format", type=str, required=False, default="short")
    parser.add_argument("-sr", '--signalp_result', type=str, required=False, default="cdhit_result")

    parser.add_argument("-t", "--tmhmm", type=bool, required=False, default=False)
    parser.add_argument("-ti", "--tmhmm_input", type=str, required=False, default="None")
    parser.add_argument("-tg", "--tmhmm_getorf", type=bool, required=False, default=False)

    args = parser.parse_args()


    args.getorf_result = path + '/' + args.getorf_result
    args.cdhit_result = path + '/' + args.cdhit_result
    args.dia_result = path + '/' + args.dia_result
    args.blastp_result = path + '/' + args.blastp_result
    filter_signalp_result =  path + '/' + 'filtered_' + args.signalp_result
    args.signalp_result = path + '/' + args.signalp_result

    start_getorf = time.time()
    if args.getorf == True:
        my_getorf(args.getorf_in, args.getorf_result, args.getorf_table, args.getorf_minsize, args.getorf_maxsize)
    else:
        print("************************* No getorf *************************")
    done_getorf = time.time()
    elapsed = done_getorf - start_getorf
    print("getorf run: " + str(elapsed/60.0))

    start_cdhit = time.time()
    if args.cdhit == True:
        if args.getorf == True:
            my_cdhit(args.getorf_result, args.cdhit_result, args.cdhit_n, args.cdhit_p, args.cdhit_c, args.cdhit_d, args.cdhit_M, args.cdhit_l, args.cdhit_s, args.cdhit_aL, args.cdhit_g)
        elif args.cdhit_input != "None":
            my_cdhit(args.cdhit_input, args.cdhit_result, args.cdhit_n, args.cdhit_p, args.cdhit_c, args.cdhit_d, args.cdhit_M, args.cdhit_l, args.cdhit_s, args.cdhit_aL, args.cdhit_g)
        else:
            print('------------------------ please input the path of the cdhit input------------------------')
    else:
        print("************************* No cdhit *************************")
    done_cdhit = time.time()
    elapsed = done_cdhit - start_cdhit
    print("cdhit run: " + str(elapsed/60.0))

    if args.diamondp == True:
        my_diamond_blastp(args.dia_db, args.cdhit_result, args.dia_result)
    else:
        print("************************* No diamond blastp *************************")

    if args.blastp == True:
        my_blastp(args.blastp_db, args.cdhit_result, args.blastp_outfmt, args.blastp_evalue, args.max_target_seqs, args.blastp_num_threads, args.blastp_result)
    else:
        print("************************* No blastp *************************")



    if args.signalp == True:
        if args.signalp_getorf == True:
            my_signalp_by_getorf(args.getorf_result, args.signalp_org, args.signalp_format, args.signalp_result)
        else:
            my_signalp_by_cdhit(args.cdhit_result, args.signalp_org, args.signalp_format, args.signalp_result)
    else:
        print("************************* No signalp *************************")



    start_tmhmm = time.time()
    if args.tmhmm == True:
        if args.tmhmm_getorf == True:
            my_tmhmm_getorf(args.getorf_result)
        else:
            if args.tmhmm_input == None:
                my_tmhmm_cdhit(path, reference=args.cdhit_result, extract_in_file_name=filter_signalp_result)
            else:
                my_tmhmm_getorf(args.tmhmm_input)
    else:
        print("************************* No tmhmm  *************************")
    done_tmhmm = time.time()
    elapsed = done_tmhmm - start_tmhmm
    print("tmhmm run: " + str(elapsed/60.0))

class myThread(threading.Thread):
    def __init__(self, threadID, name, task_name, getorf_result, direction_path):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.task_name = task_name
        self.getorf_result = getorf_result
        self.direction_path = direction_path

    def run(self):
        if self.task_name == 'tmhmm':
            start_tmhmm = time.time()
            my_tmhmm_getorf(self.getorf_result, self.direction_path)
            done_tmhmm = time.time()
            elapsed = done_tmhmm - start_tmhmm
            print("tmhmm run: " + str(elapsed / 60.0))

        elif self.task_name == 'signalp':
            start_signalp = time.time()
            my_signalp_by_getorf(self.getorf_result, signalp_org='gram-', signalp_format='short', signalp_result= self.direction_path + "/signalp_result", direction_path=self.direction_path)
            done_signalp = time.time()
            elapsed = done_signalp - start_signalp
            print("signalp run: " + str(elapsed / 60.0))

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

    args.getorf_result = path + '/' + args.getorf_result

    start_getorf = time.time()
    my_getorf(args.getorf_in, args.getorf_result, args.getorf_table, args.getorf_minsize, args.getorf_maxsize)
    done_getorf = time.time()
    elapsed = done_getorf - start_getorf
    print("getorf run: " + str(elapsed/60.0))

    # thread1 = myThread(1, 'Thread-1', 'tmhmm', args.getorf_result, args.direction_path)
    # thread2 = myThread(2, 'Thread-2', 'signalp', args.getorf_result, args.direction_path)
    #
    # thread1.start()
    # thread2.start()
    #
    # threads = []
    # threads.append(thread1)
    # threads.append(thread2)
    #
    # for t in threads:
    #     t.join()
    # print("Exiting Main Thread")

    start_tmhmm = time.time()
    my_tmhmm_getorf(args.getorf_result, args.direction_path)
    done_tmhmm = time.time()
    elapsed = done_tmhmm - start_tmhmm
    print("tmhmm run: " + str(elapsed / 60.0))

    start_signalp = time.time()
    my_signalp_by_getorf(args.getorf_result, signalp_org='gram-', signalp_format='short',
                         signalp_result=args.direction_path + "/signalp_result", direction_path=args.direction_path)
    done_signalp = time.time()
    elapsed = done_signalp - start_signalp
    print("signalp run: " + str(elapsed / 60.0))

if __name__ == '__main__':
    main()
