import unittest
import numpy as np
import glob
import os
import shutil
import time
import subprocess as sp

def cleanDir(dir):
    if os.path.exists('%s/pscond.dat' % dir):
        os.remove('%s/pscond.dat' % dir)
    if os.path.exists('%s/scond.dat' % dir):
        os.remove('%s/scond.dat' % dir)
    if os.path.exists('%s/run_magic.sh' % dir):
        os.remove('%s/run_magic.sh' % dir)
    if os.path.exists('%s/run_magic_mpi.sh' % dir):
        os.remove('%s/run_magic_mpi.sh' % dir)
    for f in glob.glob('%s/*_BIS' % dir):
        os.remove(f)
    for f in glob.glob('%s/*.test' % dir):
        os.remove(f)
    if os.path.exists('%s/stdout.out' % dir):
        os.remove('%s/stdout.out' % dir)
    for f in glob.glob('%s/*.pyc' % dir):
        os.remove(f)
    if os.path.exists('%s/__pycache__' % dir):
        shutil.rmtree('%s/__pycache__' % dir)

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
        print('\nDirectory   :           %s' % self.dir)
        print('Description :           %s' % self.description)
        self.startTime = time.time()
        cleanDir(self.dir)

        os.chdir(self.dir)
        cmd = '%s %s/input_SBDF3.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_ARS222.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_ARS443.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_CNAB2.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_SBDF4.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_PC2.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_SBDF2.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_BPR353.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_MODCNAB.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_CNLF.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_LZ232.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_KC564.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_ARS343.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_CB3.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_KC564_FD.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = '%s %s/input_KC785_FD.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

        cmd = 'cat e_kin.sbdf3 e_kin.ars222 e_kin.ars443 e_kin.cnab2 e_kin.sbdf4 e_kin.pc2 e_kin.sbdf2 e_kin.bpr353 e_kin.modcnab e_kin.cnlf e_kin.lz232 e_kin.kc564 e_kin.ars343 e_kin.cb3 e_kin.kc564_fd e_kin.kc785_fd > e_kin.test'
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))

    def tearDown(self):
        # Cleaning when leaving
        os.chdir(self.startDir)
        cleanDir(self.dir)
        for f in glob.glob('%s/*.cnab2' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.sbdf2' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.sbdf3' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.sbdf4' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.ars222' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.ars443' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.bpr353' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.pc2' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.modcnab' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.cnlf' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.lz232' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.kc564' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.ars343' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.cb3' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.kc564_fd' % self.dir):
            os.remove(f)
        for f in glob.glob('%s/*.kc785_fd' % self.dir):
            os.remove(f)

        t = time.time()-self.startTime
        st = time.strftime("%M:%S", time.gmtime(t))
        print('Time used   :                            %s' % st)

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
        datRef = readStack('%s/reference.out' % self.dir)
        datTmp = readStack('%s/e_kin.test' % self.dir)
        np.testing.assert_allclose(datRef, datTmp, rtol=self.precision, atol=1e-20)
