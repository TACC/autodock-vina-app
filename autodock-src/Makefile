idev:
	idev -p development -m 120 -n 32 -N 1

pull:
	apptainer pull docker://wjallen/autodock_vina:1.2.3
test:
	MV2_ENABLE_AFFINITY=0 MV2_SMP_USE_CMA=0 ibrun -np 32 apptainer exec autodock_vina_1.2.3.sif python3 autodock.py -r protein.pdbqt --center="11,90.5,57.5" -s "22,24,28" -m vina -d rigid -ll /scratch/projects/docking/Test-set-compressed -n 500 -f 'EMPTY'

clean:
	rm -rf configs/ output/ results*
	rm merged_results.txt
	rm autodock.log results.tar.gz

