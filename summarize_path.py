import os


def main():
    path = '/work/juancui/yinchaohe/geno_fastas/'

    top_folders = os.listdir(path)
    output_summary = []
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
                forth_level_folders = os.listdir(third_level_folder_path)
                for forth_level_folder in forth_level_folders:
                    forth_level_folder_path = third_level_folder_path + forth_level_folder + '/'
                    output_summary.append(forth_level_folder_path)
    with open('path_summary.txt', 'a') as output_file:
        output_file.writelines(output_summary)
    output_file.close()

if __name__ == '__main__':
    main()
