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
    #if os.path.exists('%s/stdout.out' % dir):
        #os.remove('%s/stdout.out' % dir)
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
        dat = np.asarray(cut, dtype='Float64')
        #out.append(dat)
        out = np.append(out, dat)
    return out

class TestTruncations(unittest.TestCase):

    def __init__(self, testName, dir, execCmd='mpirun -n 8 ../tmp/magic.exe', 
                 log=False,precision=1e-8):
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
        if log:
            self.outFile = open("%s/stdout.out" % (self.dir), 'w') 
        else: 
            self.outFile = open(os.devnull,'wb')

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
        nphis = [96, 96, 128, 128, 192, 192, 256, 256, 288, 288, 320, 320,
                 384, 384, 400, 400, 512, 512, 640, 640, 768, 768, 800,
                 800, 864, 1024]
        mincs = [1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 
                 1, 4, 1, 4, 1, 4, 4, 4]
        str = "cat "
        for k, tag in enumerate(self.tags):
            cmd = "sed -i 's/tag.*/tag         ="+'"%s"'%tag+",/g' input.nml"
            sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
            cmd = "sed -i 's/n_phi_tot.*/n_phi_tot   =%i,/g' input.nml" % nphis[k]
            sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
            cmd = "sed -i 's/minc.*/minc        =%i,/g' input.nml" % mincs[k]
            sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))

            # Run MagIC
            cmd = '%s %s/input.nml' % (self.execCmd, self.dir)
            sp.call(cmd, shell=True, stdout=self.outFile, stderr=sp.STDOUT)

            # Concatenate e_kin files
            str += 'e_kin.%s ' % tag
            cmd = "cat e_kin.%s >> e_kin.test" % tag
        cmd = str+ '> e_kin.test'
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
        self.outFile.close()

    def tearDown(self):
        # Clean up
        for tag in self.tags:
            for f in glob.glob('%s/*.%s' % (self.dir, tag)):
                os.remove(f)

        # Restore initial values in the namelist
        cmd = "sed -i 's/tag.*/tag         ="+'"test96"'+",/g' input.nml"
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
        cmd = "sed -i 's/n_phi_tot.*/n_phi_tot   =96,/g' input.nml"
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
        cmd = "sed -i 's/minc.*/minc        =1,/g' input.nml"
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))

        # Cleaning when leaving
        os.chdir(self.startDir)
        cleanDir(self.dir)

        t = time.time()-self.startTime
        st = time.strftime("%M:%S", time.gmtime(t))
        print('Time used   :                            %s' % st)

        if hasattr(self, '_outcome'): # python 3.4+
            result = self.defaultTestResult()
            self._feedErrorsToResult(result, self._outcome.errors)
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
