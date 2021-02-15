import json
import time

from pafpy import PafFile


def get_brief_records(path, cazyfamilies, sequences):
    reads = {}
    with open(path) as fileobj:
        paf = PafFile(fileobj)
        for record in paf:
            sequence = record.tname.split('|')[0]
            if sequence in sequences['name_list']:
                print("the sequence exists")
            else:
                seq_info = {}
                seq_info['seq_read_count'] = 0
                seq_info['seq_length'] = record.tlen
                sequences['name_list'][sequence] = seq_info
                sequences['seq_amount'] += 1

            multi_cazyfamilies = record.tname.split('|')[1:]
            print(multi_cazyfamilies)
            for cazyfamily in multi_cazyfamilies:
                if cazyfamily in cazyfamilies['name_list']:
                    print("the cazyfamily exists")
                    seq_array = [sequence, record.tlen]
                    if seq_array not in cazyfamilies['name_list'][cazyfamily]['cazyfamily_sequence_list']:
                        cazyfamilies['name_list'][cazyfamily]['cazyfamily_sequence_list'].append(seq_array)
                else:
                    cazy_info = {}
                    cazy_info['cazyfamily_read_count'] = 0
                    seq_array = [sequence, record.tlen]
                    cazy_info['cazyfamily_sequence_list'] = [seq_array]
                    cazyfamilies['name_list'][cazyfamily] = cazy_info
                    cazyfamilies['family_amount'] += 1

            read_id = record.qname.split('/')[0]
            if read_id in reads.keys():
                print('the read_id exists in the list')
            else:
                read = {}
                read['cazy_family'] = record.tname.split('|')[1:]
                read['seq_id'] = record.tname.split('|')[0]
                reads[read_id] = read

    output_name = path.split('.paf')[0] + '.json'
    with open(output_name, 'w') as f:
        json.dump(reads, f)
    f.close()

    return reads

def cazy_family_counter(cazyfamilies, reads_R1, reads_R2):
    for read_id in reads_R1.keys():
        if read_id in reads_R2.keys():
            read_in_R1 = reads_R1[read_id]
            for read1_cazyfamily in read_in_R1['cazy_family']:
                cazyfamilies['name_list'][read1_cazyfamily]['cazyfamily_read_count'] += 0.5
        else:
            read_in_R1 = reads_R1[read_id]
            for read1_cazyfamily in read_in_R1['cazy_family']:
                cazyfamilies['name_list'][read1_cazyfamily]['cazyfamily_read_count'] += 1

    for read_id in reads_R2.keys():
        if read_id in reads_R1.keys():
            read_in_R2 = reads_R2[read_id]
            for read2_cazyfamily in read_in_R2['cazy_family']:
                cazyfamilies['name_list'][read2_cazyfamily]['cazyfamily_read_count'] += 0.5
        else:
            read_in_R2= reads_R2[read_id]
            for read2_cazyfamily in read_in_R2['cazy_family']:
                cazyfamilies['name_list'][read2_cazyfamily]['cazyfamily_read_count'] += 1

    for cazyfamily_name in cazyfamilies['name_list']:
        cazyfamily = cazyfamilies['name_list'][cazyfamily_name]
        seq_amount = len(cazyfamily['cazyfamily_sequence_list'])
        sum_seq_length = 0
        for seq in cazyfamily['cazyfamily_sequence_list']:
            sum_seq_length += seq[1]
        average_length = sum_seq_length / seq_amount
        cazyfamily['average_seq_length'] = average_length


def sequence_counter(sequences, reads_R1, reads_R2):
    for read_id in reads_R1.keys():
        if read_id in reads_R2.keys():
            read_in_R1 = reads_R1[read_id]
            read1_sequence = read_in_R1['seq_id']
            sequences['name_list'][read1_sequence]['seq_read_count'] += 0.5
        else:
            read_in_R1 = reads_R1[read_id]
            read1_sequence = read_in_R1['seq_id']
            sequences['name_list'][read1_sequence]['seq_read_count'] += 1

    for read_id in reads_R2.keys():
        if read_id in reads_R1.keys():
            read_in_R2 = reads_R2[read_id]
            read2_sequence = read_in_R2['seq_id']
            sequences['name_list'][read2_sequence]['seq_read_count'] += 0.5
        else:
            read_in_R2 = reads_R2[read_id]
            read2_sequence = read_in_R2['seq_id']
            sequences['name_list'][read2_sequence]['seq_read_count'] += 1

def total_reads_counter(reads_R1, reads_R2):
    total_reads = 0
    for read_id in reads_R1.keys():
        if read_id in reads_R2.keys():
            total_reads += 0.5
        else:
            total_reads += 1

    for read_id in reads_R2.keys():
        if read_id in reads_R1.keys():
            total_reads += 0.5
        else:
            total_reads += 1

    return total_reads

def sequence_FPKM(sequences, amount_all_reads):
    for sequence_name in sequences['name_list'].keys():
        sequence = sequences['name_list'][sequence_name]
        convert_all_reads = amount_all_reads / pow(10, 6)
        convert_length = sequence['seq_length'] / 1000
        FPKM = sequence['seq_read_count'] / (convert_all_reads * convert_length)
        sequence['FPKM'] = FPKM

    with open('bowtie2_sequqnces.json', 'w') as f:
        json.dump(sequences, f)
    f.close()

def cazyfamily_FPKM(cazyfamilies, amount_all_reads):
    for cazyfamily_name in cazyfamilies['name_list'].keys():
        cazyfamily = cazyfamilies['name_list'][cazyfamily_name]
        convert_all_reads = amount_all_reads / pow(10, 6)
        convert_length = cazyfamily["average_seq_length"] / 1000
        FPKM = cazyfamily["cazyfamily_read_count"] / (convert_length * convert_all_reads)
        cazyfamily['FRKM'] = FPKM

    with open('bowtie2_cazyfamilies.json', 'w') as f:
        json.dump(cazyfamilies, f)
    f.close()

def main():
    cazyfamilies = {}
    cazyfamilies['name_list'] = {}
    cazyfamilies['family_amount'] = 0

    sequences = {}
    sequences['name_list'] = {}
    sequences['seq_amount'] = 0

    start = time.time()
    path_R1 = "output/b.R1.paf"
    reads_R1 = get_brief_records(path_R1, cazyfamilies, sequences)

    path_R2 = "output/b.R2.paf"
    reads_R2 = get_brief_records(path_R2, cazyfamilies, sequences)

    cazy_family_counter(cazyfamilies, reads_R1, reads_R2)
    sequence_counter(sequences, reads_R1, reads_R2)


    amount_all_reads = total_reads_counter(reads_R1, reads_R2)
    print(amount_all_reads)
    sequence_FPKM(sequences, amount_all_reads)
    cazyfamily_FPKM(cazyfamilies, amount_all_reads)
    stop = time.time()
    print('run: ' + str(stop - start))

if __name__ == '__main__':
    main()
