import os
import concurrent.futures

def parallelize_command_fill(paths):
    with concurrent.futures.ProcessPoolExecutor(max_workers=42) as executor:
        results = [executor.submit(compute_clus, _) for _ in paths]
        #results = [executor.submit(alignment_kalign, _) for _ in paths]
        for f in concurrent.futures.as_completed(results):
            print(f.result)
    print("FINISH ALL")
    return True

def compute_clus(paths):
    os.system('python3 ./scripts/transcriptOrthology.py -talg {} -gtot {} -nhxt {} -tsm {} -outf {} -lowb {}'.format(paths[0], paths[1], paths[2], paths[3], paths[4], paths[5]))
    return True
    

if __name__ == '__main__':
    folder_input = './inputs/clustering'
    for tsm_value in [1, 2, 3, 4, 5, 6]:
        directory_main = './outputs/{}'.format(tsm_value)
        os.system('mkdir {}'.format(directory_main))
        for lowb_value in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            directory_main_second = './outputs/{}/{}'.format(tsm_value, lowb_value)
            os.system('mkdir {}'.format(directory_main_second))
            entries = []
            for repository_number in range(1, 3):#53
                if repository_number not in [15, 16]:
                    mapping_file = '{}/data_{}/mappings.maps'.format(folder_input, repository_number)
                    msa_file = '{}/data_{}/msa.alg'.format(folder_input, repository_number)
                    tree_file = '{}/data_{}/tree.nhx'.format(folder_input, repository_number)

                    #apply program
                    folder_out = './outputs/{}/{}/data_{}'.format(tsm_value, lowb_value, repository_number)
                    os.system('mkdir {}'.format(folder_out))
                    entries.append((msa_file, mapping_file, tree_file, tsm_value, folder_out, lowb_value))
                    #os.system('python3 ./scripts/transcriptOrthology.py -talg {} -gtot {} -nhxt {} -tsm {} -outf {} -lowb {}'.format(msa_file, mapping_file, tree_file, tsm_value, folder_out, lowb_value))

            parallelize_command_fill(entries)