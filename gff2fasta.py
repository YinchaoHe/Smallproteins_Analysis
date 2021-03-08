import gzip
import os


def folder_maker(top_folder, sec_level_folder, third_level_folder, file):
    target_folder = 'geno_fastas'
    try:
        os.mkdir(target_folder)
    except:
        pass
    try:
        top_folder_path = target_folder + '/' + top_folder
        os.mkdir(top_folder_path)
    except:
        pass
    try:
        sec_level_folder_path = top_folder_path + '/' + sec_level_folder
        os.mkdir(sec_level_folder_path)
    except:
        pass
    try:
        third_level_folder_path = sec_level_folder_path + '/' + third_level_folder
        os.mkdir(third_level_folder_path)
    except:
        pass

    try:
        file_folder = file.split('.gff.gz')[0]
        final_folder_path = third_level_folder_path + '/' + file_folder
        os.mkdir(final_folder_path)
    except:
        pass

    return final_folder_path

def gff2fa(input_file, target_folder):
    target_file = input_file.split('/')[-1].split('.gff')[0] + '.faa'
    with gzip.open(input_file, 'rb') as f:
        all_content = f.read()
        split_tag = "##FASTA".encode()
        fasta = all_content.split(split_tag)[1]
    f.close()

    with open(target_folder + '/' + target_file, 'wb') as f:
        f.write(fasta)
    f.close()


def main():
    path = '/mnt/array2/boweny/UHGG/ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/all_genomes/'
    top_folders = ['MGYG-HGUT-001', 'MGYG-HGUT-005','MGYG-HGUT-009','MGYG-HGUT-013','MGYG-HGUT-017','MGYG-HGUT-029','MGYG-HGUT-033']
    for top_folder in top_folders:
        top_folder_path = path + top_folder + '/'
        print('first_level: ' + top_folder_path)
        sec_level_folders = os.listdir(top_folder_path)
        for sec_level_folder in sec_level_folders:
            sec_level_folder_path = top_folder_path + sec_level_folder + '/'
            print(' ' * 4 + 'sec_level_folder: ' + sec_level_folder_path)
            third_level_folders = os.listdir(sec_level_folder_path)
            for third_level_folder in third_level_folders:
                third_level_folder_path = sec_level_folder_path + third_level_folder + '/'
                print(' ' * 8 + 'third_level_folder: ' + third_level_folder_path)
                files = os.listdir(third_level_folder_path)
                for file in files:
                    file_path = third_level_folder_path + file
                    print(' ' * 12 + 'file: ' + file_path)
                    target_folder = folder_maker(top_folder, sec_level_folder, third_level_folder, file)
                    try:
                        gff2fa(file_path, target_folder)
                    except:
                        with open('error_report.txt', 'a+') as f:
                            f.write(file_path + '\n')
                        f.close()

if __name__ == '__main__':
    #main()
    file_path = 'test/GUT_GENOME160731.gff.gz'
    gff2fa(file_path, 'test')
