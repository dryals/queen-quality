#!/bin/bash

cd blup

echo "starting blup cv run $2 / 5"

./blupf90+ ${1},par1 &> cv-lastrun.log

mv solutions ../data/sol-cv${2}.txt

echo "DONE"
