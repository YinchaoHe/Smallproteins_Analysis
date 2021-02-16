import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time

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

def seq_len_signalp(signalp_result, reference):
    signalp_result = signalp_result + '_summary.signalp5'
    with open(signalp_result, 'r') as file:
       input_results = file.readlines()
    file.close()

    records = list(SeqIO.parse(reference, "fasta"))

    add_length_result = [input_results[0]]
    title = input_results[1].split("\n")[0] + " " * 4 + 'sequence_length' + '\n'
    add_length_result.append(title)

    seq_dic = {}
    for record in records:
        seq_dic[record.id] = str(len(record.seq))

    for input_result in input_results[2:]:
        id = input_result.split()[0]
        if id in seq_dic.keys():
            leng_seq = input_result.split("\n")[0] + " " * 4 + seq_dic[id] + '\n'
            add_length_result.append(leng_seq)



    signalp_result = signalp_result.split('/')[1]
    with open('intermediate/add_Seqlength_'+signalp_result, 'w') as f:

        f.writelines(add_length_result)
    f.close()

def filter_signalp(signalp_result):
    signalp_result = signalp_result + '_summary.signalp5'
    signalp_result_name_blocks = signalp_result.split('/')
    path = signalp_result_name_blocks[0]

    with open('intermediate/add_Seqlength_'+signalp_result_name_blocks[1], 'r') as file:
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


def assemble_fasta(path, extract_in_file_name, out_id_file_name, fasta_file, reference):
    extract_ids(path, extract_in_file_name, out_id_file_name)
    with open(path + '/txt_processing/' + out_id_file_name, 'r') as f:
        ids = f.readlines()
    f.close()

    search_ids = []
    for id in ids:
        id = id.split('\n')[0]
        search_ids.append(id)

    records = list(SeqIO.parse(reference, "fasta"))
    n_fasta_seqs = []
    for record in records:
        if record.id in search_ids:
#           print(record)
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

def my_getorf(getorf_in, getorf_result, getorf_table, getorf_minsize, getorf_maxsize):
    getorf_command = 'getorf -sequence ' + getorf_in + ' ' +  '-outseq' + ' ' + getorf_result + ' ' + '-table' + ' ' + getorf_table
    getorf_command = getorf_command + ' ' + '-minsize' + ' ' + getorf_minsize + ' ' + '-maxsize' + ' ' + getorf_maxsize
    os.system(getorf_command)

def my_cdhit(cdhit_input, cdhit_result, cdhit_n, cdhit_p, cdhit_c, cdhit_d, cdhit_M, cdhit_l, cdhit_s, cdhit_aL, cdhit_g):
    cdhit_command = 'cd-hit -i ' + cdhit_input + ' ' + '-o' + ' ' + cdhit_result + ' ' + '-n' + ' ' + cdhit_n + ' ' + '-p' + ' ' + cdhit_p
    cdhit_command = cdhit_command +  ' ' + '-c' + ' ' + cdhit_c + ' ' + '-d' + ' ' + cdhit_d + ' ' + '-M' + ' ' + cdhit_M + ' ' + '-l' + ' ' + cdhit_l + ' ' + '-s' + ' ' + cdhit_s + ' ' + '-aL' + ' ' + cdhit_aL + ' ' + '-g' + ' ' + cdhit_g
    os.system(cdhit_command)

def my_diamond_blastp(dia_db, cdhit_result, dia_result):
    diamond_blastp_command = 'diamond blastp --db ' + dia_db + ' ' + '-q' + ' ' + cdhit_result + ' ' + '-o' + ' ' + dia_result
    os.system(diamond_blastp_command)

def my_blastp(blastp_db, cdhit_result, blastp_outfmt, blastp_evalue, blastp_max_target_seqs, blastp_num_threads, blastp_result):
    blastp_command = 'blastp -db ' + blastp_db + ' -query ' + cdhit_result + ' -outfmt \"' + blastp_outfmt + '\" -evalue ' + blastp_evalue
    blastp_command = blastp_command + ' -max_target_seqs ' + blastp_max_target_seqs + ' -num_threads ' + blastp_num_threads + ' -out ' + blastp_result
    os.system(blastp_command)

def my_signalp_by_cdhit(cdhit_result, signalp_org, signalp_format, signalp_result):
    signalp_command = 'signalp -fasta ' + cdhit_result + ' ' + '-org' + ' ' + signalp_org + ' ' + '-format' + ' ' + signalp_format + ' ' + '-prefix' + ' ' + signalp_result
    start = time.time()
    os.system(signalp_command)
    seq_len_signalp(signalp_result=signalp_result, reference=cdhit_result)
    filter_signalp(signalp_result)
    done = time.time()
    elapsed = done - start
    print("signalp run: " + str(elapsed))

def my_signalp_by_getorf(getorf_result, signalp_org, signalp_format, signalp_result):
    signalp_command = 'signalp -fasta ' + getorf_result + ' ' + '-org' + ' ' + signalp_org + ' ' + '-format' + ' ' + signalp_format + ' ' + '-prefix' + ' ' + signalp_result
    start = time.time()
    os.system(signalp_command)
    seq_len_signalp(signalp_result=signalp_result, reference=getorf_result)
    filter_signalp(signalp_result)
    done = time.time()
    elapsed = done - start
    print("signalp run: " + str(elapsed))

def my_tmhmm_cdhit(path, reference, extract_in_file_name):
    out_id = 'chosen_ids.txt'
    fasta_file = 'chosen_seq4tmhmm.faa'
    assemble_fasta(path, extract_in_file_name=extract_in_file_name, out_id_file_name=out_id, reference=reference, fasta_file=fasta_file)
    tmhmm_input = path + '/' + fasta_file
    output = 'intermediate/tmhmm_result.txt'
    tmhmm_command = './tmhmm-2.0c/bin/tmhmm ' + tmhmm_input + ' -short > ' + output
    print("tmhmm coming.......")
    try:
        os.system(tmhmm_command)
    except:
        print("tmhmm error")

def my_tmhmm_getorf(getorf_result):
    output = 'intermediate/tmhmm_result.txt'
    tmhmm_command = './tmhmm-2.0c/bin/tmhmm ' + getorf_result + ' -short > ' + output
    print("tmhmm coming.......")
    try:
        os.system(tmhmm_command)
    except:
        print("tmhmm error")

def main():
    path = 'intermediate'
    try:
        os.mkdir(path)
    except:
        pass

    parser = argparse.ArgumentParser()
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

    parser.add_argument("-d", "--diamondp", type=bool, required=False, default=False)
    parser.add_argument("-dd", "--dia_db", type=str, required=False, default='/mnt/array2/smallproteins/database/DM_database/SmProt_KnownDatabase.dmnd')
    parser.add_argument("-dr", "--dia_result", type=str, required=False, default="cdhited_diamond_KnownDatabase.txt")

    parser.add_argument("-b", '--blastp', type=bool, required=False, default=False)
    parser.add_argument("-bd", "--blastp_db", type=str, required=False, default='/mnt/array2/smallproteins/database/Kyrpides_small_proteins/Homologs_of_All_4539_Families/cluster_homologs_db')
    parser.add_argument("-bo", "--blastp_outfmt", type=str, required=False, default='6 std qcovs qcovhsp')
    parser.add_argument("-be", '--blastp_evalue', type=str, required=False, default='10')
    parser.add_argument("-bm", "--max_target_seqs", type=str, required=False, default='1')
    parser.add_argument("-bt", "--blastp_num_threads", type=str, required=False, default='15')
    parser.add_argument("-br", "--blastp_result", type=str, required=False, default='DHIT_query_KypridesHomologs_db_evalue-le10.blastpOUT')

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
    print("getorf run: " + str(elapsed))

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
    print("cdhit run: " + str(elapsed))

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
    print("tmhmm run: " + str(elapsed))

if __name__ == '__main__':
    main()
