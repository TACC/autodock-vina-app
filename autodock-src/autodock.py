from vina import Vina
from mpi4py import MPI
import subprocess
import pickle # For unpickling pickled ligand files
import os
import blosc # For decompressing compressed ligand files
from os.path import exists
from os.path import basename # Used in the sorting function
import argparse # To accept user inputs as command line arguments
import time
import logging
import shutil


"""
Setup base MPI declarations
In this MPI implementation, rank 0 acts as the director that all other ranks 
work for. Rank 0 first completes all pre-processing. As a final step, it 
creates a list containing full filepaths to all ligand files in the ligand 
library. Then it accepts messages from worker ranks saying they are ready for 
more work. It sends work for as long as there are ligand files left. Then rank 
0 waits to hear from all ranks that they have finished processing, then 
proceeds to do the post-processing work. 
"""
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

# Initialize logging
format_str=f'[%(asctime)s {RANK}] %(filename)s:%(funcName)s:%(lineno)s - %(levelname)s: %(message)s'
logging.basicConfig(level=logging.INFO,
                    format= format_str)
                    #filename='autodock.log',
                    #filemode='w')

# Assign user inputs
parser = argparse.ArgumentParser(description ='vina commands')
parser.add_argument('-r','--receptor',type=str,
                    help='Receptor for processing')
parser.add_argument('-c','--center',type=str,required=True,
                    help='center position')
parser.add_argument('-s','--size',type=str,required=True,
                    help='box size')
parser.add_argument('-m','--module',type=str,required=True,
                    help='vina or ad4 module')
parser.add_argument('-d','--docking',type=str,required=True,
                    help='basic or flexible docking')
parser.add_argument('-f','--sidechains',type=str,default='EMPTY',required=False,
                    help='sidechains if flexible docking')
parser.add_argument('-n','--number',type=int,required=True,
                    help='top n results to return')
parser.add_argument('-ll','--ligand_library',type=str,required=True,
                    help='ligand library')
args = parser.parse_args()

# Define global constants
# Assign user inputs to grid center x-y-z and box size boundaries x-y-z
CENTER_X = float((args.center).split(",")[0])
CENTER_Y = float((args.center).split(",")[1])
CENTER_Z = float((args.center).split(",")[2])
SIZE_X = float((args.size).split(",")[0])
SIZE_Y = float((args.size).split(",")[1])
SIZE_Z = float((args.size).split(",")[2])
USER_CONFIGS = {'center_x': CENTER_X,
                'center_y': CENTER_Y,
                'center_z': CENTER_Z,
                'size_x': SIZE_X,
                'size_y': SIZE_Y,
                'size_z': SIZE_Z}

FULL_RECEPTOR=args.receptor
RECEPTOR=FULL_RECEPTOR.split('.')[0]
FLEX_RECEPTOR=f'{RECEPTOR}_flex'
DOCKING = args.docking
if DOCKING == 'rigid':
    FLEXIBLE = False
elif DOCKING == 'flexible':
    FLEXIBLE = True

 
SIDECHAINS = (args.sidechains).split('_')

#this is just a comment and another
LIBRARY_SHORT = args.ligand_library.split('/')[4]
NUMBER_OF_OUTPUTS = args.number if args.number <= 1000 else 1000

logging.debug(f"CENTER = {args.center}; BOX = {args.size}; DOCKING = {DOCKING}")
logging.debug(f"LIGAND LIBRARY = {args.ligand_library}; LIBRARY SHORT = {LIBRARY_SHORT}")

# Internal constants
# tasks should be nodes * 128 / cpus
# These values were determined through internal benchmarks to allow an entire set to run
# under 24 hours
TASKS = int(os.environ['SLURM_NTASKS']) # What the user chose on the web portal
NODES = int(os.environ['SLURM_NNODES']) # What the user chose on the web portal
if LIBRARY_SHORT in ['Test-set-compressed', 'Enamine-PC-compressed', 'ZINC-fragments-compressed', 'ZINC-in-trials-compressed']:
    EXPECTED_NODES = 1
    EXPECTED_TASKS = 32
elif LIBRARY_SHORT == 'Enamine-HTSC-compressed':
    EXPECTED_NODES = 10
    EXPECTED_TASKS = 320
elif LIBRARY_SHORT == 'Enamine-AC-compressed':
    EXPECTED_NODES = 3
    EXPECTED_TASKS = 96
CPUS = 4
VERBOSITY = 0 # Prints vina docking progress to stdout if set to 1 (normal) or 2 (verbose)
POSES = 1 # If set to 1, only saves the best pose/score to the output ligand .pdbqt file
EXHAUSTIVENESS = 8
MAX_SIDECHAINS = 6

logging.debug(f"TASKS = {TASKS}; NODES = {NODES}; EXPECTED TASKS = {EXPECTED_TASKS}; EXPECTED NODES = {EXPECTED_NODES}")



def main():
    current_responses = 0
    if RANK == 0:
        start_time = time.time()
        # Pre-Processing
        pre_processing()            #step 1
        ligands = prep_ligands()    #step 2
        # Let other ranks know pre-processing is finished; they can now ask for work
        logging.info("Pre-processing finished; Starting sendrecv with all ranks")
        
        for i in range(1,SIZE):
            #sending the pre-processing message to other messages 
            COMM.send('pre-processing finished; ask for work', dest=i)
        logging.info("Send finished; All ranks sent;")
        
        pre_responses = 0
        while pre_responses != (SIZE - 1):
            response = COMM.recv(source = MPI.ANY_SOURCE)
            if response == 'ready to go':
                pre_responses += 1
        logging.info("recv finished; All ranks responded; Starting main work")

        COMM.send("Rank 1 will be ready to clean as we go", dest=1, tag=1)
        logging.info("rank 1 will be ready to clean as we go")

        # Until all ligands have been docked, send more work to worker ranks
        while ligands:
            source = COMM.recv(source=MPI.ANY_SOURCE)
            COMM.send(ligands.pop(), dest=source)
        logging.info("Rank 0: List of ligands is now empty")

        # When all ligands have been sent, let worker ranks know they can stop
        logging.info("Tell all ranks there is no more work")
        for i in range(2,SIZE):
                COMM.send('no more ligands', dest=i)
                logging.info('no more ligands')
        
        current_responses = 0
        while current_responses != (SIZE - 2):
            response = COMM.recv(source = MPI.ANY_SOURCE)
            if response == 'message received--proceed to post-processing':
                current_responses += 1

        current_time = time.strftime("%H:%M:%S")
        logging.info(f"Rank {RANK}: All ranks have responded; Proceeding to post-processing at {current_time}")
        print(f"From Rank 0: All ranks done; going to post-processing at {time.time()}")
        
        COMM.send('stop working', dest=1)
        finished_message = COMM.recv(source=1,tag=1)
        logging.info(finished_message)
        top_ligand_filenames =  COMM.recv(source=1,tag=2)

        # Post-Processing
        sort()
        isolate_output_for_rank1(top_ligand_filenames)
        reset()
        end_time = time.time()
        total_time = end_time - start_time

        logging.info(f"Script runtime = {total_time}.")
        logging.info(f"Nodes: {NODES}.")
        logging.info(f"Tasks: {TASKS}.")
        logging.info(f"Library: {LIBRARY_SHORT}.")

    elif RANK == 1:
        COMM.recv(source=0) # Wait for rank 0 to finish pre-processing
        COMM.send('ready to go', dest=0)
        message = COMM.recv(source = 0)
       
        if message == "Rank 1 will be ready to clean as we go":
            logging.info("This is means it works")
            clean_as_we_go_another()

    else: # All ranks besides rank 0
        COMM.recv(source=0) # Wait for rank 0 to finish pre-processing
        COMM.send('ready to go', dest=0)
        processing()    

# End def main()

def clean_as_we_go_another():
    '''
        cleaning as we go function constantly checks for result files being generated, sorts it using sort() function
        and keeps only the top n results by finding a threshold this gets updated 
    '''
    logging.info("We are cleaning as we go")
    ligands_folder = 'output/results/ligands'
    top_results = args.number
    file_counter = 0
    files_before_cleanup = 10
    while True:
        # receive message from any source
        message = COMM.recv(source=MPI.ANY_SOURCE)
        # Check if it's a stop message from rank 0
        if message == 'stop working':
            logging.info("Rank 0 sent a message to rank 1 to stop")
            command = "find output/pdbqt/ -name '*pdbqt' | wc -l"
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
            pdbqt_files = result.stdout.strip()
            pdbqt_files = int(pdbqt_files)
            logging.info(f'number of files in system is: {pdbqt_files}')
            #checks if we are keeping the top n elements before break
            if pdbqt_files == top_results:
                logging.info("Stopping cleaning as we go, ready for last sort")
                break
            else:
                # if they are not the same we delete the remaining files not in the top n
                logging.info("number of pdbqt_files is not the same as top_results")
                sort_for_rank1()
                flines_to_remove,lines_to_keep = file_deletion(top_results)
                command = "find output/pdbqt/ -name '*pdbqt' | wc -l"
                result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
                pdbqt_files = result.stdout.strip()
                pdbqt_files = int(pdbqt_files)
                logging.info(f'number of files in system is: {pdbqt_files}')
                break
        # check if its a file result update
        elif "File results" in message:
            file_counter+=1
            current_time = time.strftime("%H:%M:%S")
            sender_rank = message.split('_')[1]  
            logging.info(message + " at " + current_time)
            if file_counter >= files_before_cleanup:
                sort_for_rank1() # this cats the results that were generated into "merged_results" then sorts them into "sorted_scores_all.txt". 
                lines_to_remove,lines_to_keep = file_deletion(top_results)
                file_counter = 0   
    
    #informs rank 0 that rank 1 is done 
    COMM.send("finished--proceed to post-processing",dest = 0,tag=1)
    COMM.send(lines_to_keep,dest = 0,tag=2)
    logging.info("Rank 1 finished cleaning and sent completion message to Rank 0")

def file_deletion(top_results):
    with open('./output/results/sorted_scores_all.txt', 'r') as sorted_scores:
        lines = sorted_scores.readlines()
        logging.info(f'number of lines = {str(len(lines))}')
            
    # Remove Top N elements
    lines_to_remove = lines[top_results:]
    lines_to_keep = lines[:top_results]
    logging.info(f'number of lines in linestoremove= {str(len(lines_to_remove))}') #checks size for sorted_lines which should only have until the nth element
                
    try:
            for dirpath, _, filenames in os.walk('./output/pdbqt'):
                    for line in lines_to_remove:
                        name = line.split()[0]
                        for filename in filenames:
                            file_path = os.path.join(dirpath, filename)
                            if filename == f'output_{name}': # this is correct
                                os.remove(file_path)
                                logging.info(f'file path {file_path} was removed')
    except:
            logging.info("file couldnt be deleted")
                
    return lines_to_remove,lines_to_keep  
    
def check_user_configs():
    # User inputted box size must be within bounds specified below
    for size in [SIZE_X, SIZE_Y, SIZE_Z]:
        if not (size <= 30 and size >= 1):
            logging.error("box size is outside the bounds (1-30).")
            COMM.Abort()
    # User must input a file ending in .pdb or .pdbqt
    if not (FULL_RECEPTOR.endswith('.pdb') or FULL_RECEPTOR.endswith('.pdbqt')):
        logging.error("Please provide a .pdb or .pdbqt file.")
        COMM.Abort()
          
    # User inputted grid center must be within the receptor's min/max bounds
    # User inputted sidechains must be found within the .pdbqt receptor file
    # The coordinates are parsed according to the documentation on pdb files
    # found here:
    # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    all_sidechains = []
    xbounds = []
    ybounds = []
    zbounds = []
    with open(f'{RECEPTOR}.pdbqt', 'r') as r:
        line = r.readline()
        while line:
            line = r.readline()
            if line.startswith('ATOM') or line.startswith('HETATM'):
                xbounds.append(float(line[30:38]))
                ybounds.append(float(line[38:46]))
                zbounds.append(float(line[46:54]))
                all_sidechains.append(line[17:20] + line[22:26])
    if FLEXIBLE == True:
        for sidechain in SIDECHAINS:
            if not sidechain in all_sidechains:
                logging.error("Please provide valid flexible sidechain \
                              names, separated by underscores (e.g. THR315_GLU268).")
                COMM.Abort()

    logging.debug(f'X-bounds = {min(xbounds)} - {max(xbounds)} | Y-bounds = {min(ybounds)} - {max(ybounds)} | Z-bounds = {min(zbounds)} - {max(zbounds)}')       
    if not (min(xbounds)) <= CENTER_X <= (max(xbounds)):
        logging.error("Center x coordinate is not within bounds.")
        COMM.Abort()
    if not (min(ybounds)) <= CENTER_Y <= (max(ybounds)):
        logging.error("Center y coordinate is not within bounds.")
        COMM.Abort()
    if not (min(zbounds)) <= CENTER_Z <= (max(zbounds)):
        logging.error("Center z coordinate is not within bounds.")
        COMM.Abort()
    
    
    # Maximum of 6 sidechains
    if len(SIDECHAINS) > MAX_SIDECHAINS:
        logging.error("Too many sidechains specified (max: 6).")
        COMM.Abort()
    
    # User inputted #Nodes and #Tasks must match our internal values (specified above) exactly
    if not (TASKS == EXPECTED_TASKS) or not (NODES == EXPECTED_NODES):
        logging.error(f'Incorrect values for #Nodes and/or #ProcessorsPerNode.\n \
                        Please review input guidelines before submitting a job.\n \
                        Current #Nodes={NODES}.\n \
                        Current#Tasks={TASKS}.\n \
                        Expected #Nodes for {LIBRARY_SHORT}={EXPECTED_NODES}.\n \
                        Expected #Tasks for {LIBRARY_SHORT}={EXPECTED_TASKS}.')
        COMM.Abort()
    logging.debug("User configs checked. No aborts called")


def prep_config():
    # Write user inputted configurations (grid center and box size) to a 
    #   config file for AutoDock Vina to use
    with open('./configs/config.config', 'w+') as f:
        for config, value in USER_CONFIGS.items():
            f.write(f'{config} = {value}\n')
    logging.debug("Config file written")


def prep_maps():
    # Generate affinity maps using write-gpf.py and AFDRSuite's autogrid4. 
    #   Generates maps for all possible ligand atom types for a given receptor
    if args.module == 'ad4':
        if exists(f'{RECEPTOR}.gpf'):
            os.remove(f'{RECEPTOR}.gpf')
        subprocess.run([f"python3 /autodock-src/scripts/write-gpf.py --box {'./configs/config.config'} \
                        {RECEPTOR}.pdbqt"], shell=True)
        subprocess.run([f"autogrid4 -p {RECEPTOR}.gpf"], shell=True)
        logging.debug("Affinity maps generated")


def prep_receptor():
    # Converts a PDB receptor to a PDBQT, if needed. If the user has specified 
    #   flexible docking, also prepares the rigid receptor and user-chosen 
    #   flexible sidechains.
    if FULL_RECEPTOR.endswith('.pdb'):
        try:
            subprocess.run([f'prepare_receptor -r {RECEPTOR}.pdb'], shell=True)
            logging.debug("Receptor prepped successfully")
        except Exception as e:
            logging.error(f"error on rank {RANK}: error prepping receptor")
            logging.debug(e)
            COMM.Abort()
    if FLEXIBLE == True:
        try:
            subprocess.run([f"pythonsh /autodock-src/scripts/prepare_flexreceptor.py \
                -g {RECEPTOR}.pdbqt -r {RECEPTOR}.pdbqt \
                -s {'_'.join(SIDECHAINS)}"], shell=True)
            logging.debug("Flex receptor prepped successfully")
        except Exception as e:
            logging.error(f"error on rank {RANK}: error prepping flex receptor")
            logging.debug(e)
            COMM.Abort()
    

def prep_ligands():
    """
        Function prep_ligands returns a list where each item is the path to a pickled and compressed text file 
        containing multiple ligand strings (ignores pkl or .dat format).

        Return: 
            ligand_paths (list): a list where each item containis multiple ligand strings (the items arepickled and compressed).
    
    """
    
    ligand_paths = []
    for dirpath, _, filenames in os.walk(args.ligand_library):
        for filename in filenames:
            if filename.endswith('.pkl') or filename.endswith('.dat'):
                ligand_paths.append(f'{dirpath}/{filename}')
    logging.info(f"Ligands prepped. Number of ligand batch files = {len(ligand_paths)}")
    return ligand_paths


def run_docking(ligands, v, directory): #step 3 
    # Runs AutoDock on each ligand in the given set; outputs a .pdbqt file 
    #   showing the pose and all scores; appends the ligand name (filename) 
    #   and it's best pose/score to a temporary results file
    if not bool(ligands):
        logging.error("Ligand batch was set to empty; no ligands docked")
        return
    output_directory = f'./output/pdbqt/{RANK}{directory}'
    if not exists(output_directory):
        os.makedirs(output_directory)
    for _, filename in enumerate(ligands):
        ligand = ligands[filename]
        if not bool(ligand) or ligand == {}:
            logging.error(f"Rank {RANK}: Ligand {filename} string was empty; skipping")
            continue
        try:
            v.set_ligand_from_string(ligand)
        except Exception as e:
            logging.error(f"Rank {RANK}: Error setting ligand {filename}; Error = {e}")
            continue
        try:
            v.dock(exhaustiveness=EXHAUSTIVENESS)
        except Exception as e:
            logging.error(f"Rank {RANK}: Error docking ligand {filename}; Error = {e}")
            continue
        try:
            v.write_poses(f'{output_directory}/output_{filename}', \
                           n_poses=POSES, overwrite=True)
        except Exception as e:
            logging.error(f"Rank {RANK}: Error writing ligand {filename}; Error = {e}")
            continue
        ### putting subprocess write function on one line to avoid race condition
        subprocess.run([f"grep -i -m 1 'REMARK VINA RESULT:' \
                       {output_directory}/output_{filename} \
                       | awk '{{print $4}}' >> results_{RANK}.txt; echo {filename} \
                       >> results_{RANK}.txt"], shell=True)
        COMM.send(f"File results_{RANK}.txt was generated", dest=1)

def unpickle_and_decompress(path_to_file):
    # Given a filepath, decompresses and unpickles the file. Returns the 
    #   contents of the file as a dictionary where keys are ligand filenames 
    #   and values are the actual ligand strings
    with open(path_to_file, 'rb') as f:
        compressed_pickle = f.read()
    try:
        depressed_pickle = blosc.decompress(compressed_pickle)
    except Exception as e:
            logging.error(f"Error on rank {RANK}: Error decompressing ligand batch {compressed_pickle}; Error = {e}")
            depressed_pickle = ''
    try:
        dictionary_of_ligands = pickle.loads(depressed_pickle)
    except Exception as e:
            logging.error(f"Error on rank {RANK}: Error unpickling ligand batch {depressed_pickle}; Error = {e}")
            dictionary_of_ligands = {}
    return dictionary_of_ligands


def pre_processing():
    """
        pre_processing functions that helps to reduce clutter in main(), 
        it prepares the config files, the receptor, affinity maps and checks for user configuration
    
    """
    
    # Helper function to reduce clutter in main()
    os.makedirs('configs')
    os.makedirs('output/pdbqt')
    os.makedirs('output/results/ligands')
    prep_config()
    prep_receptor()
    prep_maps()
    check_user_configs()


def processing():
    """
        Function processing handles sending ligands to worker ranks, receiving ligands from rank 0, and running docking for each ligand set.
    """
    # Long method; difficult to de-couple any one part without breaking many things

    # Initialize docking configurations
    # Note: Internal benchmarks showed that on a given node with 128 cores, 
    #   setting ibrun -np to 32 and, below, setting CPU=4 granted the fastest 
    #   docking. Verbosity set to 0 to increase speeds, as stdout cannot be 
    #   captured from Vina's Python module.
    if args.module == 'vina':
        v = Vina(sf_name='vina', cpu=CPUS, verbosity=VERBOSITY)
        if FLEXIBLE == True:
            v.set_receptor(f'{RECEPTOR}.pdbqt', f'{FLEX_RECEPTOR}.pdbqt')
        else:
            v.set_receptor(f'{RECEPTOR}.pdbqt')
        v.compute_vina_maps(center=[CENTER_X, CENTER_Y, CENTER_Z],
                            box_size=[SIZE_X, SIZE_Y, SIZE_Z])
    elif args.module == 'ad4':
        v = Vina(sf_name='ad4', cpu=CPUS, verbosity=VERBOSITY)
        v.load_maps(map_prefix_filename = RECEPTOR)
        
    # Ask rank 0 for ligands and dock until rank 0 says done
    count = 1
    directory = 1
    while True:
        COMM.send(RANK,dest=0) # Ask rank 0 for another set of ligands ###
        ligand_set_path = COMM.recv(source=0) # Wait for a response 

       # COMM.send(RANK,dest = 1)

        if ligand_set_path == 'no more ligands':
            current_time = time.strftime("%H:%M:%S")
            logging.info(f"Rank {RANK} has finished all work at {current_time}")
            print(f"Rank {RANK} done at {time.time()}")
            COMM.send('message received--proceed to post-processing',dest=0)
            break
        try:
            ligands = unpickle_and_decompress(ligand_set_path)
        except Exception as e:
            logging.error(f'Error on rank {RANK}: could not \
                            unpickle/decompress {ligand_set_path}.')
            logging.debug(e)
        try:
            run_docking(ligands, v, directory)
        except Exception as e:
            #logging.error(f'Error on rank {RANK}: docking error with \
            #                ligand set {ligand_set_path}, ligands {ligands}.')
            logging.error(f'Error on rank {RANK}: docking error with ligand set {ligand_set_path}.')
            logging.debug(e)

        count += 1
        if count == 100:
            count = 1
            directory += 1


def sort(): # step 4 
    """
        The sort() function cats all results files into one, it arranges ligands based on the highest score; prints these sorted results are written to 
        sorted_scores.txt; finally cleans up the directory
    """
    # Cats all results files into one, arranges each line to read: 
    #   (ligand, top score), then sorts by score so that highest scoring 
    #   ligands are on top; prints these sorted results are written to 
    #   sorted_scores.txt; finally cleans up the directory

    subprocess.run(["cat results* >> merged_results.txt"], shell=True)
    INPUTFILE = 'merged_results.txt'
    OUTPUTFILE = './output/results/sorted_scores.txt'
    result = []

    with open(INPUTFILE) as data:
        line = data.readline()
        while line:
            filename = basename(line.split()[-1])
            v = data.readline().split()[0]
            result.append(f'{v} {filename}\n')
            line = data.readline()

    with open(OUTPUTFILE, 'w') as data:
        data.writelines(sorted(result[:NUMBER_OF_OUTPUTS], \
                        key=lambda x: float(x.split()[1])))
        
def sort_for_rank1(): # step 4 
    """
        The sort() function cats all results files into one, it arranges ligands based on the highest score; prints these sorted results are written to 
        sorted_scores.txt; finally cleans up the directory
    """
    # Cats all results files into one, arranges each line to read: 
    #   (ligand, top score), then sorts by score so that highest scoring 
    #   ligands are on top; prints these sorted results are written to 
    #   sorted_scores.txt; finally cleans up the directory

    try:
        os.remove('./output/results/sorted_scores_all.txt')
        os.remove('merged_results.txt')
    except:
        logging.info('exception in sort_for_rank1')
    subprocess.run(["cat results* >> merged_results.txt"], shell=True)
    INPUTFILE = 'merged_results.txt'
    OUTPUTFILE = './output/results/sorted_scores_all.txt'
    result = []

    with open(INPUTFILE) as data:
        line = data.readline()    
        while line:
            filename = basename(line.split()[-1])
            #logging.info(filename)
            v = data.readline().split()[0]
            #logging.info(v)
            result.append(f'{v} {filename}\n')
            line = data.readline()

    with open(OUTPUTFILE, 'w') as data:
        data.writelines(sorted(result, key=lambda x: float(x.split()[1])))

def isolate_output_for_rank1(top_ligand_filenames):
   
    logging.info("runing isolate_output_for_1")
    top_ligand_filenames = [line.split()[0] for line in top_ligand_filenames]
    ligands_docked = 0
    
    for dirpath, _, filenames in os.walk('./output/pdbqt'):
        ligands_docked += len(filenames)
        for top_filename in top_ligand_filenames:
            for filename in filenames:
                if filename == f'output_{top_filename}':
                    shutil.move(f'{dirpath}/{filename}', './output/results/ligands')
    
    combined = open('./output/results/combined_docked_ligands.pdbqt', 'w+')
    while top_ligand_filenames:
        with open(f'./output/results/ligands/output_{top_ligand_filenames.pop(0)}') as f:
            lines = f.read()
            combined.write(lines)
    combined.close()
    
    logging.info(f"Total ligands docked = {ligands_docked}")
    subprocess.run([f'tar -czf results.tar.gz ./output/results'], shell=True)

def isolate_output(): #step 5 isolates
    """
        isolate_ouput() function copies the user-specified top n ligand output files to a single directory (top_ligand_filenames)
    """
    # Copies the user-specified top n ligand output files to a single directory
    top_ligand_filenames = []
    ligands_docked = 0
    logging.info("runing isolate_output")
    with open('./output/results/sorted_scores.txt', 'r') as results:
        for _, line in enumerate(results):
            top_ligand_filenames.append(line.split()[0])

    for dirpath, _, filenames in os.walk('./output/pdbqt'):
        ligands_docked += len(filenames)
        for top_filename in top_ligand_filenames:
            for filename in filenames:
                if filename == f'output_{top_filename}':
                    shutil.move(f'{dirpath}/{filename}', './output/results/ligands')
    
    combined = open('./output/results/combined_docked_ligands.pdbqt', 'w+')
    while top_ligand_filenames:
        with open(f'./output/results/ligands/output_{top_ligand_filenames.pop(0)}') as f:
            lines = f.read()
            combined.write(lines)
    combined.close()
    
    logging.info(f"Total ligands docked = {ligands_docked}")
    subprocess.run([f'tar -czf results.tar.gz ./output/results'], shell=True)


def reset():
    """
     reset() function resets directory by removing all created files from this job
    
    """

    for dirpath, _, filenames in os.walk('.'):
        for filename in filenames:
            if filename.endswith(('.map', '.txt', '.gpf', '.fld', '.xyz')):
                logging.debug(f"Removed {filename}\n")
                os.remove(f'{dirpath}/{filename}') 
    shutil.rmtree('./output')
    shutil.rmtree('./configs')
        
main()
