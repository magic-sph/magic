#!/usr/bin/env python
import numpy as np
import os
import re
import getpass
import string
import sys
from subprocess import check_output
'''
This script is used to resubmit # jobs with dependencies in the current working
directory. This only works with SBATCH systems.

./resub.py 10
'''

uname = getpass.getuser()
pwd = os.getcwd()


if len(sys.argv) > 1:
    n_resub = int(sys.argv[1])
else:
    n_resub = 1

def get_old_tag(pbs_file):
    tag = check_output("grep -oP 'tag\s+=\s?+\K(.*),' %s" % pbs_file, shell=True)
    tag = tag.decode('utf8').strip('\n')
    tag = tag.split('"')[1] # Get rid of the ""

    return tag


def submit_job(dep, jobfile):
    job_sub = check_output("sbatch --dependency=afterok:%i %s" % (dep, jobfile),
                           shell=True)
    job_sub = job_sub.decode('utf8').strip('\n')
    job_sub = int(job_sub.split()[-1])

    return job_sub


def get_job_dir(jobname):
    run_path = check_output("scontrol show jobid %i | grep -oP 'WorkDir=\K(.*)'" % jobname,
                            shell=True)
    run_path = run_path.decode('utf8').strip('\n')
    run_path = os.path.realpath(run_path)

    if run_path == pwd:
        path = True
    else:
        path = False

    return path


def get_job_file(jobname, path):

    cmd = check_output("scontrol show jobid %i | grep -oP 'Command=%s\K(.*)'" % (jobname, path),
                            shell=True)
    cmd = cmd.decode('utf8').strip('\n')[1:]

    return cmd

integers = re.compile('[0-9]')
batchFiles = re.compile(r'run([0-9]+)\.sbatch')
tagName = re.compile(r' *tag *= *"([\w\.]*)" *, *')
    
cmd = "squeue -u %s | awk '{print $1}' > list_jobs" % uname
os.system(cmd)
jobs = np.loadtxt('list_jobs', skiprows=1)
jobs = np.sort(jobs)

# Determine the jobs that are already in the queue
job_cwd = np.array([], dtype=np.int16)
for job in jobs:

    if get_job_dir(int(job)):
        job_cwd = np.append(job_cwd, int(job))

# Determine last job name
most_recent = job_cwd[-1]
old_job = get_job_file(most_recent, pwd)

for i in range(n_resub):

    # Determine new job file
    if batchFiles.match(old_job):
        oldindex = int(batchFiles.match(old_job).group(1))
        new_job = 'run%i.sbatch' % (oldindex+1)
    
    # Determine the new tag
    old_tag = get_old_tag(old_job)
    if integers.match(old_tag[-1]) or old_tag[-1] == 'Z':
        new_tag = oldtag + 'a'
    else:
        ind = string.ascii_letters.index(old_tag[-1])
        new_tag = old_tag[:-1] + string.ascii_letters[ind+1]
    
    # Determine the new checkpoint file
    new_checkpoint = 'checkpoint_end.%s' % old_tag
    
    # Now create the job file
    cmd = 'cp %s %s' % (old_job, new_job)
    os.system(cmd)
    
    # SED substitutions
    cmd = "sed -i 's/ *tag *=.*/tag         =\"%s\",/g' %s" % (new_tag, new_job)
    os.system(cmd)
    cmd = "sed -i 's/ *start_file.*/start_file  =\"%s\",/g' %s" % (new_checkpoint, new_job)
    os.system(cmd)
    cmd = "sed -i 's/ *l_start_file.*/l_start_file=.true.,/g' %s" % new_job
    os.system(cmd)
    cmd = "sed -i 's/ *init_u.*/init_u     =0,/g' %s" % new_job
    os.system(cmd)
    cmd = "sed -i 's/ *init_t.*/init_t     =0,/g' %s" % new_job
    os.system(cmd)
    
    # Now submit the new job and overwrite old by new
    most_recent = submit_job(most_recent, new_job)
    old_job = new_job
    print(old_job)


# Remove useless file
os.system('rm list_jobs')
