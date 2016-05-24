#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe openmpi-smp 8
#$ -q *@@3rd_gen

source /etc/profile.d/valet.sh

vpkg_require "vasp/5.3.2+vtst+d3-pgi10"

echo "GridEngine parameters:"
echo "  nhosts         = $NHOSTS"
echo "  nproc          = $NSLOTS"
echo "  mpiexec        =" `which mpiexec`
echo "  pe_hostfile    = $PE_HOSTFILE"
echo "  vasp           =" `which vasp`
echo
cat $PE_HOSTFILE
echo
echo "-- begin OPENMPI run --"
time mpiexec --n $NSLOTS --host localhost --mca btl sm,self --display-map vasp
echo "-- end OPENMPI run --"
