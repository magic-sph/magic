import unittest
import numpy as np
import glob
import os
import shutil
import time
import subprocess as sp

def cleanDir(dir):
    if os.path.exists('{}/pscond.dat'.format(dir)):
        os.remove('{}/pscond.dat'.format(dir))
    if os.path.exists('{}/scond.dat'.format(dir)):
        os.remove('{}/scond.dat'.format(dir))
    if os.path.exists('{}/run_magic.sh'.format(dir)):
        os.remove('{}/run_magic.sh'.format(dir))
    if os.path.exists('{}/run_magic_mpi.sh'.format(dir)):
        os.remove('{}/run_magic_mpi.sh'.format(dir))
    for f in glob.glob('{}/*_BIS'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.test'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.cnab2'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.sbdf2'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.sbdf3'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.sbdf4'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.ars222'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.ars443'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.bpr353'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.pc2'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.modcnab'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.cnlf'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.lz232'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.kc564'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.ars343'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.cb3'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.kc564_fd'.format(dir)):
        os.remove(f)
    for f in glob.glob('{}/*.kc785_fd'.format(dir)):
        os.remove(f)
    if os.path.exists('{}/stdout.out'.format(dir)):
        os.remove('{}/stdout.out'.format(dir))
    for f in glob.glob('{}/*.pyc'.format(dir)):
        os.remove(f)
    if os.path.exists('{}/__pycache__'.format(dir)):
        shutil.rmtree('{}/__pycache__'.format(dir))

def readStack(file):
    f = open(file, 'r')
    #out = []
    out = np.array([])
    for line in f.readlines():
        cut = line.split()
        dat = np.asarray(cut, dtype=np.float64)
        #out.append(dat)
        out = np.append(out, dat)
    return out

class TimeSchemes(unittest.TestCase):

    def __init__(self, testName, dir, execCmd='mpirun -n 8 ../tmp/magic.exe', 
                 precision=1e-8):
        super(TimeSchemes, self).__init__(testName)
        self.dir = dir
        self.precision = precision
        self.execCmd = execCmd
        self.startDir = os.getcwd()
        self.description = "Test several time integrators"

    def list2reason(self, exc_list):
        if exc_list and exc_list[-1][0] is self:
            return exc_list[-1][1]

    def setUp(self):
        # Cleaning when entering
        print('\nDirectory   :           {}'.format(self.dir))
        print('Description :           {}'.format(self.description))
        self.startTime = time.time()
        cleanDir(self.dir)

        os.chdir(self.dir)
        cmd = '{} {}/input_SBDF3.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_ARS222.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_ARS443.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_CNAB2.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_SBDF4.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_PC2.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_SBDF2.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_BPR353.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_MODCNAB.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_CNLF.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_LZ232.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_KC564.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_ARS343.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_CB3.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_KC564_FD.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '{} {}/input_KC785_FD.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = 'cat e_kin.sbdf3 e_kin.ars222 e_kin.ars443 e_kin.cnab2 e_kin.sbdf4 e_kin.pc2 e_kin.sbdf2 e_kin.bpr353 e_kin.modcnab e_kin.cnlf e_kin.lz232 e_kin.kc564 e_kin.ars343 e_kin.cb3 e_kin.kc564_fd e_kin.kc785_fd > e_kin.test'
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))

    def tearDown(self):
        # Cleaning when leaving
        os.chdir(self.startDir)
        cleanDir(self.dir)

        t = time.time()-self.startTime
        st = time.strftime("%M:%S", time.gmtime(t))
        print('Time used   :                            {}'.format(st))

        if hasattr(self, '_outcome'): # python 3.4+
            if hasattr(self._outcome, 'errors'):  # python 3.4-3.10
                result = self.defaultTestResult()
                self._feedErrorsToResult(result, self._outcome.errors)
            else:  # python 3.11+
                result = self._outcome.result
        else:  # python 2.7-3.3
            result = getattr(self, '_outcomeForDoCleanups', 
                             self._resultForDoCleanups)

        error = self.list2reason(result.errors)
        failure = self.list2reason(result.failures)
        ok = not error and not failure

        if ok:
            print('Validating results..                     OK')
        else:
            if error:
                print('Validating results..                     ERROR!')
                print('\n')
                print(result.errors[-1][-1])
            if failure:
                print('Validating results..                     FAIL!')
                print('\n')
                print(result.failures[-1][-1])

    def outputFileDiff(self):
        datRef = readStack('{}/reference.out'.format(self.dir))
        datTmp = readStack('{}/e_kin.test'.format(self.dir))
        np.testing.assert_allclose(datRef, datTmp, rtol=self.precision, atol=1e-20)
