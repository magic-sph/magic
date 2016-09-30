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
        dat = np.asarray(cut, dtype='Float64')
        #out.append(dat)
        out = np.append(out, dat)
    return out


def generateEkinFile(fileName='e_kin.test'):
    from magic import MagicGraph, Movie

    gr = MagicGraph(ivar=1)

    file = open(fileName, 'w')
    st = '%.4f %.4f %.4f %.4f %.4f %.4f %.4f' % ( gr.entropy[10, 13, 3], 
         gr.Br[33, 0, 2], gr.Btheta[3, 11, 11], gr.Bphi[34, 12, 25], 
         gr.vr[11, 15, 2], gr.vtheta[33, 33, 3], gr.vphi[32, 33, 7] )

    # Write output for graphic files
    file.write(st+'\n')

    av = Movie(file='AV_mov.start', iplot=False)
    ahf = Movie(file='AHF_mov.start', iplot=False)
    brcmb = Movie(file='Br_CMB_mov.start', iplot=False)
    vtr = Movie(file='Vt_R=C1_mov.start', iplot=False)
    vreq = Movie(file='Vr_EQU_mov.start', iplot=False)
    teq = Movie(file='T_EQU_mov.start', iplot=False)
    vortz = Movie(file='VorZ_EQU_mov.start', iplot=False)
    hel = Movie(file='HE_R=C2_mov.start', iplot=False)
    st = '%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f' % (av.data[1, 121, 12], 
         ahf.data[0, 99, 33], brcmb.data[1, 47, 128], vtr.data[0, 33, 87],
         vreq.data[0, 28, 31], teq.data[1, 29, 29], vortz.data[1, 1, 1],
         hel.data[0, 3, 4])

    # Write output for movie files
    file.write(st)
    file.close()


class TestGraphicMovieOutputs(unittest.TestCase):

    def __init__(self, testName, dir, execCmd='mpirun -n 8 ../tmp/magic.exe'):
        super(TestGraphicMovieOutputs, self).__init__(testName)
        self.dir = dir
        self.execCmd = execCmd
        self.startDir = os.getcwd()
        self.description = "Test Graphic and Movie outputs"

    def list2reason(self, exc_list):
        if exc_list and exc_list[-1][0] is self:
            return exc_list[-1][1]

    def setUp(self):
        # Cleaning when entering
        print('\nDirectory   :           %s' % self.dir)
        print('Description :           %s' % self.description)
        self.startTime = time.time()
        cleanDir(self.dir)
        for f in glob.glob('%s/*.start' % self.dir):
            os.remove(f)
        os.chdir(self.dir)
        cmd = '%s %s/input.nml' % (self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

    def tearDown(self):
        # Cleaning when leaving
        os.chdir(self.startDir)
        cleanDir(self.dir)
        for f in glob.glob('%s/*.start' % self.dir):
            os.remove(f)

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

    @unittest.skipUnless('MAGIC_HOME' in os.environ, 
                         'MAGIC_HOME is not defined! source sourceme.sh!')
    def outputFileDiff(self):
        generateEkinFile('e_kin.test')
        datRef = readStack('%s/reference.out' % self.dir)
        datTmp = readStack('%s/e_kin.test' % self.dir)
        np.testing.assert_equal(datRef, datTmp)
