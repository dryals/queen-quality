#!/bin/bash
ulimit -s unlimited

echo "starting blup cv run $2 / 5"

./blupf90+ blup/${1} &> blup/cv-lastrun.log

mv solutions sol-cv${2}.txt

echo "DONE"
