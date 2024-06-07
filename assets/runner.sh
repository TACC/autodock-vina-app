#!/bin/bash

# set flexible sidechains, if any
FLEX=${flexible_sidechains}
if [ -z $FLEX ]; then
    FLEX='Empty'
fi

# mpi settings - set VAL based on number of nodes
if [[ "${_tapisNodes}" > 1 ]]; then
    VAL=1
else 
    VAL=0
fi

# Log commands, timing, run job
echo -n "starting: "
date

echo "================================================================"
echo "python3 /autodock-src/autodock.py \\"
echo "    -r ${receptor} \\"
echo "    --center=\"${center_x},${center_y},${center_z}\" \\"
echo "    -s \"${size_x},${size_y},${size_z}\" \\"
echo "    -m ${forcefield} \\"
echo "    -d ${docking} \\"
echo "    -ll ${library} \\"
echo "    -n ${top_n_scores} \\"
echo "    -f ${FLEX} "
echo "================================================================"

MV2_ENABLE_AFFINITY=0 MV2_SMP_USE_CMA=${VAL} \
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
