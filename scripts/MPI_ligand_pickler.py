"""This is an alternative script for pickling and compressing ligands that utilizes MPI for parallel processing. It runs significantly faster than the previous script, but the output directories and filenames are less clearly organized (though still working). A script to rename/reorganize the output is in order."""

import os
import os.path
import pickle
import blosc
import copy
from mpi4py import MPI
import subprocess

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

MAX_LIGANDS_PER_SET = 100
path = '/scratch/02875/docking/test/Enamine-PC-test'
write_path = '/scratch/02875/docking/test/Enamine-HTSC-compressed2'
write_path = './test'
ligands = {}

def main():
    if rank == 0:
        names = []
        for dirpath, dirnames, filenames in os.walk(path):
            for filename in [f for f in filenames]:    
                names.append(f'{dirpath}/{filename}')
        for i in range(1,size):
            comm.sendrecv('pre-processing finished; ask for work', dest=i)
        while names:
            lig_set=[]
            for i in range(100):
                if names:
                    lig_set.append(names.pop())
            source = comm.recv(source=MPI.ANY_SOURCE)
            comm.send(lig_set, dest=source)
        for i in range(1,size):
            comm.send('done', dest=i)
            comm.recv(source=i)
    else:
        comm.recv(source=0)
        comm.send(rank, dest=0)
        count = 0
        while True:
            comm.send(rank, dest=0)
            lig_set = comm.recv(source=0)
            if lig_set == 'done':
                comm.send('ok', dest=0)
                break
            ligands = {}
            for lig in lig_set:
                ligand = open(lig, 'r')
                filename = lig.split('/')[-1]
                ligands[filename] = ligand.read()
            pickled_dict = pickle.dumps(ligands)
            compressed_pickle = blosc.compress(pickled_dict)
            if not os.path.exists(f'{write_path}/{rank}'):
                os.makedirs(f'{write_path}/{rank}')
            with open(f'{write_path}/{rank}/ligands_{rank}_{count}.pkl', 'wb') as f:
                f.write(compressed_pickle)
            count += 1


main()
