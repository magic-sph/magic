#!/bin/sh

echo "Verify against reference.out"
python3 compare.py

echo "Clean dir"

rm -r magic.exe* slurm-* *.dat *.test* Magic_*
