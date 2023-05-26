from __future__ import print_function
import unittest
import numpy as np
import glob
import os
import time
import shutil
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
    for f in glob.glob('%s/*.start' % dir):
        os.remove(f)
    for f in glob.glob('%s/*.continue' % dir):
        os.remove(f)
    if os.path.exists('%s/stdout.out' % dir):
        os.remove('%s/stdout.out' % dir)
    for f in glob.glob('%s/*.pyc' % dir):
        os.remove(f)
    if os.path.exists('%s/__pycache__' % dir):
        shutil.rmtree('%s/__pycache__' % dir)


def readData(file):
    return np.loadtxt(file)


class ConductingRotatingIC(unittest.TestCase):

    def __init__(self, testName, dir, execCmd='mpirun -n 8 ../tmp/magic.exe', 
                 precision=1e-8):
        super(ConductingRotatingIC, self).__init__(testName)
        self.dir = dir
        self.precision = precision
        self.execCmd = execCmd
        self.startDir = os.getcwd()
        self.description = "Conducting and rotating inner core"

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
        cmd = '%s %s/input.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))
        cmd = '%s %s/input_restart.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))
        cmd = 'cat e_kin.start e_kin.continue > e_kin.test'
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))
        cmd = 'cat e_mag_oc.start e_mag_oc.continue > e_mag_oc.test'
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

    def tearDown(self):
        # Cleaning when leaving
        os.chdir(self.startDir)
        cleanDir(self.dir)

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
        # Kinetic energy
        datRef = readData('%s/reference.out' % self.dir)
        datTmp = readData('%s/e_kin.test' % self.dir)
        np.testing.assert_allclose(datRef, datTmp, rtol=self.precision, 
                                   atol=1e-20)

        # Magnetic energy
        datRef = readData('%s/referenceMag.out' % self.dir)
        datTmp = readData('%s/e_mag_oc.test' % self.dir)
        np.testing.assert_allclose(datRef, datTmp, rtol=self.precision,
                                   atol=1e-20)
