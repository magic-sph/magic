#!/bin/sh

echo "Verify results"

python3 compare.py

echo "Clean dir"

rm -r magic.exe* slurm-* *.dat *.test*
