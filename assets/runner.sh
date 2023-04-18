 Allow over-ride
if [ -z "${CONTAINER_IMAGE}" ]
then
    version=$(cat ./_util/VERSION)
    CONTAINER_IMAGE="index.docker.io/library/ubuntu:bionic"
fi
. lib/container_exec.sh

# Log commands, timing, run job
echo -n "starting: "
date

echo "================================================================"
echo "singularity pull vina_1.2.3.0.sif docker://austindarrow/autodock_vina:1.2.3.0"
echo "================================================================"
echo "MV2_ENABLE_AFFINITY=0 ibrun singularity exec vina_1.2.3.0.sif python3 autodock.py"
echo "-r ${receptor}"
echo "--center="${center_x},${center_y},${center_z}""
echo "-s "${size_x},${size_y},${size_z}""
echo "-m ${forcefield}"
echo "-d ${docking}"
echo "-ll ${library}"
echo "-n ${top_n_scores}"
echo "-f ${flex}"
echo "================================================================"

singularity pull vina_1.2.3.0.sif docker://austindarrow/autodock_vina:1.2.3.0

flex=${flexible_sidechains}
if
     [ -z $flex ]
then
     flex='Empty'
fi
    
MV2_ENABLE_AFFINITY=0 ibrun singularity exec vina_1.2.3.0.sif python3 autodock.py \
     -r ${receptor} \
     --center="${center_x},${center_y},${center_z}" \
     -s "${size_x},${size_y},${size_z}" \
     -m ${forcefield} \
     -d ${docking} \
     -ll ${library} \
     -n ${top_n_scores} \
     -f ${flex}

echo -n "ending: "
date