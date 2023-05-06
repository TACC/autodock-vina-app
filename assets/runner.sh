# Allow over-ride
if [ -z "${CONTAINER_IMAGE}" ]
then
    CONTAINER_IMAGE="wjallen/autodock_vina:1.2.3"
fi

# Silence xalt errors
module unload xalt
apptainer pull --disable-cache vina_1.2.3.sif docker://${CONTAINER_IMAGE}


# set flexible sidechains, if any
FLEX=${flexible_sidechains}
if
    [ -z $FLEX ]
then
    FLEX='Empty'
fi

# mpi settings
MPI_SETTING=""
if [ "${library}" = "/scratch/projects/docking/Test-set-compressed" ]; then         
    MPI_SETTING="0"               
elif [ "${library}" = "/scratch/projects/docking/Enamine-PC-compressed" ]; then     
    MPI_SETTING="0"               
elif [ "${library}" = "/scratch/projects/docking/Enamine-AC-compressed" ]; then     
    MPI_SETTING="1"                                 
elif [ "${library}" = "/scratch/projects/docking/Enamine-HTSC-compressed" ]; then   
    MPI_SETTING="1"                                 
elif [ "${library}" = "/scratch/projects/docking/ZINC-fragments-compressed" ]; then 
    MPI_SETTING="0"               
elif [ "${library}" = "/scratch/projects/docking/ZINC-in-trials-compressed" ]; then 
    MPI_SETTING="0"               
fi      


# Log commands, timing, run job
echo -n "starting: "
date

echo "================================================================"
echo "MV2_ENABLE_AFFINITY=0 MV2_SMP_USE_CMA=${MPI_SETTING} ibrun apptainer exec vina_1.2.3.sif \ "
echo "    python3 /autodock-src/autodock.py \ "
echo "    -r ${receptor} \ "
echo "    --center="${center_x},${center_y},${center_z}" \ "
echo "    -s "${size_x},${size_y},${size_z}" \ "
echo "    -m ${forcefield} \ "
echo "    -d ${docking} \ "
echo "    -ll ${library} \ "
echo "    -n ${top_n_scores} \ "
echo "    -f ${FLEX} "
echo "================================================================"

MV2_ENABLE_AFFINITY=0 MV2_SMP_USE_CMA=${MPI_SETTING} ibrun apptainer exec vina_1.2.3.sif \
     python3 /autodock-src/autodock.py \
     -r ${receptor} \
     --center="${center_x},${center_y},${center_z}" \
     -s "${size_x},${size_y},${size_z}" \
     -m ${forcefield} \
     -d ${docking} \
     -ll ${library} \
     -n ${top_n_scores} \
     -f ${FLEX}

echo -n "ending: "
date

# clean up
rm vina_1.2.3.sif

