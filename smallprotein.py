import argparse
import os
import assem_fasta

def my_getorf(getorf_in, getorf_result, getorf_table, getorf_minsize, getorf_maxsize):
    getorf_command = 'getorf -sequence ' + getorf_in + ' ' +  '-outseq' + ' ' + getorf_result + ' ' + '-table' + ' ' + getorf_table
    getorf_command = getorf_command + ' ' + '-minsize' + ' ' + getorf_minsize + ' ' + '-maxsize' + ' ' + getorf_maxsize
    os.system(getorf_command)

def my_cdhit(getorf_result, cdhit_result, cdhit_n, cdhit_p, cdhit_c, cdhit_d, cdhit_M, cdhit_l, cdhit_s, cdhit_aL, cdhit_g):
    cdhit_command = 'cd-hit -i ' + getorf_result + ' ' + '-o' + ' ' + cdhit_result + ' ' + '-n' + ' ' + cdhit_n + ' ' + '-p' + ' ' + cdhit_p
    cdhit_command = cdhit_command +  ' ' + '-c' + ' ' + cdhit_c + ' ' + '-d' + ' ' + cdhit_d + ' ' + '-M' + ' ' + cdhit_M + ' ' + '-l' + ' ' + cdhit_l + ' ' + '-s' + ' ' + cdhit_s + ' ' + '-aL' + ' ' + cdhit_aL + ' ' + '-g' + ' ' + cdhit_g
    os.system(cdhit_command)

def my_diamond_blastp(dia_db, cdhit_result, dia_result):
    diamond_blastp_command = 'diamond blastp --db ' + dia_db + ' ' + '-q' + ' ' + cdhit_result + ' ' + '-o' + ' ' + dia_result
    os.system(diamond_blastp_command)

def my_signalp(cdhit_result, signalp_org, signalp_format, signalp_result):
    signalp_command = 'signalp -fasta ' + cdhit_result + ' ' + '-org' + ' ' + signalp_org + ' ' + '-format' + ' ' + signalp_format + ' ' + '-prefix' + ' ' + signalp_result
    os.system(signalp_command)
    path = signalp_result.split('/')[0]
    signalp_result = signalp_result + '_summary.signalp5'
    assem_fasta.filter_signalp(path, signalp_result)

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
    parser.add_argument("-db", "--dia_db", type=str, required=False, default='/mnt/array2/smallproteins/database/DM_database/SmProt_KnownDatabase.dmnd ')
    parser.add_argument("-dr", "--dia_result", type=str, required=False, default="cdhited_diamond_KnownDatabase.txt")

    parser.add_argument("-s", "--signalp", type=bool, required=False, default=False)
    parser.add_argument("-so", "--signalp_org", type=str, required=False, default="gram-" )
    parser.add_argument("-sf", "--signalp_format", type=str, required=False, default="short")
    parser.add_argument("-sr", '--signalp_result', type=str, required=False, default="cdhit_result")
    args = parser.parse_args()


    args.getorf_result = path + '/' + args.getorf_result
    args.cdhit_result = path + '/' + args.cdhit_result
    args.dia_result = path + '/' + args.dia_result
    args.signalp_result = path + '/' + args.signalp_result

    if args.getorf == True:
        my_getorf(args.getorf_in, args.getorf_result, args.getorf_table, args.getorf_minsize, args.getorf_maxsize)
    else:
        print("************************* No getorf *************************")

    if args.cdhit == True:
        my_cdhit(args.getorf_result, args.cdhit_result, args.cdhit_n, args.cdhit_p, args.cdhit_c, args.cdhit_d, args.cdhit_M, args.cdhit_l, args.cdhit_s, args.cdhit_aL, args.cdhit_g)
    else:
        print("************************* No cdhit *************************")

    if args.diamondp == True:
        my_diamond_blastp(args.dia_db, args.cdhit_result, args.dia_result)
    else:
        print("************************* No diamond blastp *************************")

    if args.signalp == True:
        my_signalp(args.cdhit_result, args.signalp_org, args.signalp_format, args.signalp_result)
    else:
        print("************************* No signalp *************************")




if __name__ == '__main__':
    main()
