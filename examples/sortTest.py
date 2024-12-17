from mpi4py import MPI
import time
import random
import os
import math
import subprocess

# MPI comm declarations
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

"""
This script is designed to mimic the control flow of autodock.py
to run on your machine:
mpiexec -n {number_of_processes} python sortTest.py

We're comparing two different versions of autodock.py:
    sort_then_trim() mimics the current implementation on Lonestar6
    remember_MVS() is an altered version designed to speed up score sorting
"""


def prep_ligands():
    """
    instead of passing actual ligands and docking them, we're
    pretending! 1000 elements in ligand path, each with 100
    random integers
    """
    ligand_paths = []
    for i in range(1000):
        ligand = []
        for i in range(100):
            element = random.randint(1, 9999)
            ligand.append(element)
        ligand_paths.append(ligand)

    print("ligands prepped!")

    return ligand_paths


def sort_for_rank1():
    print("R1 is sorting...")
    try:
        os.remove("./output/results/sorted_scores_all.txt")
    except:
        pass
    try:
        os.remove("merged_results.txt")
    except:
        pass

    subprocess.run(["touch merged_results.txt"], shell=True)
    subprocess.run(["cat results_*.txt >> merged_results.txt"], shell=True)
    inputfile = "merged_results.txt"
    outputfile = "./output/results/sorted_scores_all.txt"
    results = []

    with open(inputfile, "r") as fin:
        line = fin.readline()
        while line:
            if line != "\n":
                results.append(f"{line}")
            line = fin.readline()

    subprocess.run([f"touch {outputfile}"], shell=True)
    with open(outputfile, "w") as fout:
        fout.writelines(sorted(results, key=lambda x: float(x)))
    print(f"R1 sorted {len(results)} results")

    return


def delete_files(num):
    with open("./output/results/sorted_scores_all.txt", "r") as fin:
        lines = fin.readlines()

    lines_to_remove = lines[num:]
    lines_to_keep = lines[:num]

    print(f"Rank 1 deleted {len(lines_to_remove)} files after cleaning")
    return lines_to_keep


def clean_as_we_go():
    """
    here, we read a bunch of results from the workers, sort them
    periodically, and keep the top 1000 scores. The detritus gets dumped.
    """
    file_counter = 0
    lines_to_keep = []

    while True:
        # receive message from any source
        message = COMM.recv(source=MPI.ANY_SOURCE)

        # check if R0 told me to stop
        if message == "STOP":
            print("R0 told R1 to stop cleaning")
            break
        elif "File results" in message:
            file_counter += 1
            if file_counter >= 1000:
                print(f"file counter hit {file_counter}, starting clean up")
                sort_for_rank1()
                lines_to_keep = delete_files(1000)
                file_counter = 0

    # tell R0 we're done cleaning
    COMM.send("stopping clean_as_we_go", dest=0, tag=1)
    COMM.send(lines_to_keep, dest=0, tag=2)

    return


def run_docking(ligands, directory):
    """
    pretend we're docking. we wait a small random amount of time to simulate
    the docking process. we do NOT sort the data as it comes, we write it
    in the order we receive it
    """
    if not bool(ligands):
        print("ligand batch was set to empty: no ligands docked")
        return
    # in the real version. we write each score to the following directory
    output_directory = f"./output/pdbqt/{RANK}{directory}"
    # if not exists(output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    # and then merge all those directories into one results_{RANK}.txt
    # instead we're saving the 'scores' to a list and writing the list
    # directly to results_{RANK}.txt

    scores = []
    for ligand in ligands:
        # v.dock(exhaustiveness = EXHAUSTIVENESS)
        docking_score = math.sqrt(ligand)
        # time.sleep((random.randint(1, 30)) / 10)
        # v.write_poses(f'{output_directory}/output_{filename})
        """
        rather than write each score to its own file and merge those
        scores to results_{RANK}.txt, we'll add each score to an array
        and write that to disk at the end
        """
        scores.append(docking_score)
    with open(f"results_{RANK}.txt", "a") as fout:
        fout.write("\n".join(str(i) for i in scores))
        fout.write("\n")

    print(f"R{RANK} wrote {len(scores)} to results_{RANK}.txt")
    COMM.send(f"File results_{RANK} was generated", dest=1)


def processing():
    """
    a facsimile of the vina docking job. we'll find the square root
    of the number we received, save it to an array, and wait a small
    random amount of time
    """
    # ask rank 0 for ligands and 'dock' until rank 0 says done
    print(f"R{RANK} entered processing loop")
    count = 1
    directory = 1
    while True:
        COMM.send(RANK, dest=0)  # tell root where to send ligand data
        ligands = COMM.recv(source=0)
        # print(f"R{RANK} received {ligands}")

        if ligands == "no more ligands":
            print(f"R{RANK} has no more work to do")
            COMM.send("proceed to post-processing", dest=0)
            break

        try:
            print(f"R{RANK} is docking...")
            run_docking(ligands, directory)
        except Exception as e:
            print(f"ERROR: {e}")

        count += 1
        if count == 100:
            count = 1
            directory += 1


def sort_then_trim():
    """
    RANK 0: controller/root
    RANK 1: clean as we go
    RANK 2-N: pretend work (dock molecules in the real guy)
    """
    if RANK == 0:
        start_time = time.time()
        # pretend we're pre-processing
        ligands = prep_ligands()

        # let other ranks know pre-processing is finished; they ask for work
        print("instructing workers to ask for work")
        for i in range(1, SIZE):
            COMM.send("request work", dest=i)

        pre_responses = 0
        while pre_responses != (SIZE - 1):
            response = COMM.recv(source=MPI.ANY_SOURCE)
            if response == "ready":
                print("R0 received ready ACK")
                pre_responses += 1

        print("rank 0 received ACK from workers, telling R1 to clean")
        COMM.send("Rank 1 will be ready to clean as we go", dest=1, tag=1)

        # dispatch each ligand to worker ranks
        while ligands:
            print("R0 dispatching ligands...")
            source = COMM.recv(source=MPI.ANY_SOURCE)
            print(f"R0 heard ACK from R{source}")
            payload = ligands[0]
            ligands = ligands[1:]
            COMM.send(payload, dest=source)

        # when all ligands are sent, let worker ranks know they can stop
        print("no more ligands")
        for i in range(2, SIZE):
            COMM.send("no more ligands", dest=i)

        current_responses = 0
        while current_responses != (SIZE - 2):
            response = COMM.recv(source=MPI.ANY_SOURCE)
            if response == "proceed to post-processing":
                current_responses += 1

        print("rank 0 received ACK stop from workers, telling R1 to stop working")
        # tell rank 1 to stop working
        COMM.send("STOP", dest=1)
        finished_message = COMM.recv(source=1, tag=1)
        top_ligands_message = COMM.recv(source=1, tag=2)
        # print(top_ligands_message)
        end_time = time.time()
        total_time = end_time - start_time
        print(
            f"sort_then_trim finished; job was finished in {'\033[92m'}{total_time}{'\033[0m'} seconds"
        )
        print(f"Rank 0 received {len(top_ligands_message)} with the first five values:")
        for i in range(5):
            print(top_ligands_message[i])

    elif RANK == 1:
        COMM.recv(source=0)  # wait for rank 0 to finish pre-processing
        print("rank 1 heard that R0 is done pre-processing, sending ready message")
        COMM.send("ready", dest=0)
        message = COMM.recv(source=0)
        if message == "Rank 1 will be ready to clean as we go":
            print("rank 1 is cleaning")
            clean_as_we_go()

    else:
        COMM.recv(source=0)  # wait for rank 0 to finish pre-processing
        print(f"R{RANK} heard that R0 is finished pre-processing")
        COMM.send("ready", dest=0)  # worker rank is ready to receive input data
        processing()

    return


# --- DIFFERENT APPROACH: remember_MVS() ---


def sort_with_MVS(top_scores, scores):
    # print("R1 is sorting...")
    top_scores += scores
    results_sorted = sorted(top_scores, key=lambda x: float(x))

    return results_sorted[:1000]


def clean_and_update_MVS(minimum_viable_score):
    lines_to_keep = []

    while True:
        message = COMM.recv(source=MPI.ANY_SOURCE)

        if message == "STOP":
            print("R0 told R1 to stop cleaning")
            break
        elif "File results" in message[0]:
            scores = message[1]
            if len(scores) > 0:
                lines_to_keep = sort_with_MVS(lines_to_keep, scores)
                # lines_to_keep = delete_files(1000)

            minimum_viable_score = lines_to_keep[-1]
            # print(f"sending new MVS ({minimum_viable_score}) to R0")
            COMM.send(minimum_viable_score, dest=0)

    print("R1 broke out of cleaning loop")
    # tell R0 we're done cleaning
    COMM.send("stopping clean_and_update_MVS", dest=0, tag=1)
    COMM.send(lines_to_keep, dest=0, tag=2)

    return


def run_docking_MVS(ligands, minimum_viable_score):
    if not bool(ligands):
        return

    scores = []
    for ligand in ligands:
        # v.dock(exhaustiveness = EXHAUSTIVENESS)
        docking_score = math.sqrt(ligand)
        # v.write_poses(f'{output_directory}/output_{filename}')
        if docking_score < float(minimum_viable_score):
            scores.append(docking_score)

    with open(f"results_{RANK}.txt", "a") as fout:
        fout.write("\n".join(str(i) for i in scores))
        fout.write("\n")

    print(f"R{RANK} wrote {len(scores)} to results_{RANK}.txt")
    payload = [f"File results_{RANK} was generated", scores]
    COMM.send(payload, dest=1)
    print(f"told R1 that R{RANK} finished working")
    return


def processing_MVS():
    """
    a facsimile of the vina docking job. we find the square root of
    the int we received, and compare it with the minimum viable score
    sent from R0. if the "score" (square root) is smaller than the
    minimum viable score, we save it to an array. once all ints have
    been evaluated, write the trimmed result to disk
    """
    while True:
        COMM.send(RANK, dest=0)  # tell root where to send ligand data
        response = COMM.recv(source=0)

        if response == "no more ligands":
            print(f"R{RANK} has no more work to do")
            COMM.send("proceed to post-processing", dest=0)
            break

        minimum_viable_score = response[0]
        ligands = response[1]
        run_docking_MVS(ligands, minimum_viable_score)


def remember_MVS():
    """
    RANK 0: controller/root
    RANK 1: clean as we go
    RANK 2-N: pretend work (dock molecules in the real app)
    """
    if RANK == 0:
        start_time = time.time()
        # pretend we're pre-processing
        subprocess.run(["rm results_*.txt"], shell=True)
        ligands = prep_ligands()

        # let other ranks know pre-processing is finished; they ask for work
        print("instructiong workers to ask for work")
        for i in range(1, SIZE):
            COMM.send("request work", dest=i)

        pre_responses = 0
        while pre_responses != (SIZE - 1):
            response = COMM.recv(source=MPI.ANY_SOURCE)
            if response == "ready":
                print("R0 received ready message")
                pre_responses += 1

        print("R0 received ACK from workers, telling R1 to clean")
        COMM.send("Rank 1 will be ready to clean as we go", dest=1, tag=1)

        # dispatch each ligand to worker ranks
        # initialize a minimum viable score as a basis to determine which docking scores
        # write to disk. this value gets updated every time rank 1 sorts results
        minimum_viable_score = 100
        while ligands:
            print(f"R0 dispatching ligands ({len(ligands)} remaining)...")
            source = COMM.recv(source=MPI.ANY_SOURCE)
            payload = [minimum_viable_score, ligands[0]]
            ligands = ligands[1:]
            COMM.send(payload, dest=source)
            minimum_viable_score = COMM.recv(source=1)
            # print("R0 received new MVS")
            # minimum_viable_score = COMM.recv(source=1)

        # when all ligands are sent, let worker ranks know they can stop
        print("no more ligands")
        for i in range(2, SIZE):
            COMM.send("no more ligands", dest=i)

        current_responses = 0
        while current_responses != (SIZE - 2):
            response = COMM.recv(source=MPI.ANY_SOURCE)
            if response == "proceed to post-processing":
                current_responses += 1

        # tell rank 1 to stop working
        COMM.send("STOP", dest=1)
        top_ligands_message = COMM.recv(source=1, tag=2)
        end_time = time.time()
        total_time = end_time - start_time
        print(
            f"remember_MVS finished; job was finished in {'\033[92m'}{total_time}{'\033[0m'} seconds"
        )
        print(f"Rank 0 received {len(top_ligands_message)} with the first five values:")
        for i in range(5):
            print(top_ligands_message[i])

    elif RANK == 1:
        COMM.recv(source=0)  # wait for rank 0 to finish pre-processing
        COMM.send("ready", dest=0)
        message = COMM.recv(source=0)
        if message == "Rank 1 will be ready to clean as we go":
            clean_and_update_MVS(100)

    else:
        COMM.recv(source=0)
        COMM.send("ready", dest=0)
        processing_MVS()


def main():
    sort_then_trim()
    time.sleep(5)
    remember_MVS()


if __name__ == "__main__":
    main()
