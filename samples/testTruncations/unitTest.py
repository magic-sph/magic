import unittest
import numpy as np
import glob
import os
import time
import shutil
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

class TestTruncations(unittest.TestCase):

    def __init__(self, testName, dir, execCmd='mpirun -n 8 ../tmp/magic.exe', 
                 precision=1e-8):
        super(TestTruncations, self).__init__(testName)
        self.dir = dir
        self.precision = precision
        self.execCmd = execCmd
        self.startDir = os.getcwd()
        self.description = "Test various truncations"
        self.tags = ['test96', 'test96m4', 'test128', 'test128m4', 'test192',
                     'test192m4', 'test256', 'test256m4', 'test288', 'test288m4',
                     'test320', 'test320m4', 'test384', 'test384m4', 'test400',
                     'test400m4', 'test512', 'test512m4', 'test640', 'test640m4',
                     'test768', 'test768m4', 'test800', 'test800m4', 'test864m4',
                     'test1024m4']

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
        nphis = [96, 96, 128, 128, 192, 192, 256, 256, 288, 288, 320, 320,
                 384, 384, 400, 400, 512, 512, 640, 640, 768, 768, 800,
                 800, 864, 1024]
        mincs = [1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 
                 1, 4, 1, 4, 1, 4, 4, 4]
        str = "cat "
        for k, tag in enumerate(self.tags):
            cmd = "sed -i 's/tag.*/tag         ="+'"{}"'.format(tag)+",/g' input.nml"
            sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
            cmd = "sed -i 's/n_phi_tot.*/n_phi_tot   ={:d},/g' input.nml".format(nphis[k])
            sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
            cmd = "sed -i 's/minc.*/minc        ={:d},/g' input.nml".format(mincs[k])
            sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))

            # Run MagIC
            cmd = '{} {}/input.nml'.format(self.execCmd, self.dir)
            sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                    stderr=open(os.devnull, 'wb'))

            # Concatenate e_kin files
            str += 'e_kin.{} '.format(tag)
            cmd = "cat e_kin.{} >> e_kin.test".format(tag)
        cmd = str+ '> e_kin.test'
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))

    def tearDown(self):
        # Clean up
        for tag in self.tags:
            for f in glob.glob('{}/*.{}'.format(self.dir, tag)):
                os.remove(f)

        # Restore initial values in the namelist
        cmd = "sed -i 's/tag.*/tag         ="+'"test288"'+",/g' input.nml"
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
        cmd = "sed -i 's/n_phi_tot.*/n_phi_tot   =288,/g' input.nml"
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
        cmd = "sed -i 's/minc.*/minc        =1,/g' input.nml"
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))

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
