#!/bin/bash

python3 "part1.py" $1 $2
python "part2.py"
python3 "part3.py" $1 $2

rm "temp1.csv" "temp2.csv" "temp3.csv"