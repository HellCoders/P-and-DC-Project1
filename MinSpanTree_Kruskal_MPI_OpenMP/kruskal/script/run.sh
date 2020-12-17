#!/bin/bash

if [ ! -f g01.bin ] || [ ! -f g05.bin ] || [ ! -f g10.bin ] || [ ! -f g20.bin ]; then
    echo "Dataset(s) missing! Check /script directory."
    exit 0
fi

echo -e "SERIAL (1%):"
sres[0]="$(../bin/kruskal-ser g01.bin | grep -o '[0-9]\.\?[0-9]*')"
echo -e "SERIAL (5%):"
sres[1]="$(../bin/kruskal-ser g05.bin | grep -o '[0-9]\.\?[0-9]*')"
echo -e "SERIAL (10%):"
sres[2]="$(../bin/kruskal-ser g10.bin | grep -o '[0-9]\.\?[0-9]*')"
echo -e "SERIAL (20%):"
sres[3]="$(../bin/kruskal-ser g20.bin | grep -o '[0-9]\.\?[0-9]*')"


n=8
echo -e "MPI (1%): $n procs"
pres[0]="$(mpirun -n $n ../bin/kruskal-mpi g01.bin | grep -o '[0-9]\.\?[0-9]*')"
echo -e "MPI (5%): $n procs"
pres[1]="$(mpirun -n $n ../bin/kruskal-mpi g05.bin | grep -o '[0-9]\.\?[0-9]*')"
echo -e "MPI (10%): $n procs"
pres[2]="$(mpirun -n $n ../bin/kruskal-mpi g10.bin | grep -o '[0-9]\.\?[0-9]*')"
echo -e "MPI (20%): $n procs"
pres[3]="$(mpirun -n $n ../bin/kruskal-mpi g20.bin | grep -o '[0-9]\.\?[0-9]*')"

for i in `seq 0 3`;
do
    echo
    sarr=(${sres[i]}) # MST length, total time, time no I/O.
    parr=(${pres[i]})
    if [ "${sarr[0]}"="${parr[0]}" ];
    then
        printf "Test %d: Correct length (%d)\n" $i ${parr[0]}
        printf "Speedup       : %.2f\n" $(bc -l<<< ${sarr[1]}/${parr[1]})
        printf "Speedup no I/O: %.2f\n" $(bc -l<<< ${sarr[2]}/${parr[2]})
    fi
done
