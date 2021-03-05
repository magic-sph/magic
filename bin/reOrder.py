#!/usr/bin/env python
from magic import MagicSetup
import glob
import os

'''
This script is used to reorder a directory in which the output files from
MagIC have lost their original writing time. It will sort the log files and
then touch all the output files to recreate a consistent series.
'''

logFiles = sorted(glob.glob('log.*'))

# Determine the last log
lastLog = logFiles[-1]
restarted_from = []
for log in logFiles:
    stp = MagicSetup(nml=log, quiet=True)
    if stp.l_start_file == 'T':
        l_start_file = True
    else:
        l_start_file = False
    if l_start_file:
        tag = stp.tag
        lst = stp.start_file.split('.')
        oldtag = ''.join(lst[1:])
        if os.path.exists('log.'+oldtag): # all logFile exists
            restarted_from.append('log.'+oldtag)
diff = list(set(logFiles)-set(restarted_from))
if len(diff) == 1:
    lastLog = diff[0]
else:
    x = raw_input('Enter last tag: ')
    lastLog = 'log.' + x

print('Latest log file found is {}'.format(lastLog))

# Sort log file by reading checkpoint name to get the previous tag
logs_sorted = [lastLog]
while 1:
    stp = MagicSetup(nml=logs_sorted[0], quiet=True)
    if stp.l_start_file == 'T':
        l_start_file = True
    else:
        l_start_file = False
    if l_start_file:
        tag = stp.tag
        lst = stp.start_file.split('.')
        oldtag = ''.join(lst[1:])

        if os.path.exists('log.'+oldtag): # all logFile exists
            logs_sorted.insert(0, 'log.'+oldtag)
        else:
            break
    else:
        break

# Now touch the files with sorted tags
for log in logs_sorted:
    lst = log.split('.')
    tag = ''.join(lst[1:])
    os.system('touch *.{}'.format(tag))
