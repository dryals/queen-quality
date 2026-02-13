#!/bin/bash

cd ~/ryals/queen-quality/blup

echo "starting blup cv run $2 / 5"

./blupf90+ ${1} &> cv-lastrun.log

mv solutions ../data/sol-cv${2}.txt

echo "DONE"
