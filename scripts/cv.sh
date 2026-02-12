#!/bin/bash

cd /home/dylan/Documents/bees/harpurlab/project/gensel/data2025/blup                    

export LD_LIBRARY_PATH=~/anaconda3/lib

ulimit -s unlimited

echo "starting blup cv run $2 / 5"

./blupf90+ $1 &> cv-lastrun.log

mv solutions cv-H-${2}-sol.txt

echo "DONE"
