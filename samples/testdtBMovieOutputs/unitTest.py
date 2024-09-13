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


def generateEkinFile(fileName='e_kin.test'):
    from magic import Movie

    ab = Movie(file='AB_mov.start', iplot=False)
    at = Movie(file='ATmov.start', iplot=False)
    abadv = Movie(file='ABAdv_mov.start', iplot=False)
    abpro = Movie(file='ABPro_mov.start', iplot=False)
    brpror= Movie(file='BrPro_R=1_038461_mov.start', iplot=False)
    fladv = Movie(file='FLAdv_mov.start', iplot=False)
    fldif = Movie(file='FLDif_mov.start', iplot=False)
    flpro = Movie(file='FLPro_mov.start', iplot=False)
    jrpro3d = Movie(file='JrPro_3D_mov.start', iplot=False)
    jradvr = Movie(file='JrAdv_R=0_788461_mov.start', iplot=False)
    jrdifr = Movie(file='JrDif_R=1_288461_mov.start', iplot=False)

    file = open(fileName, 'w')
    st = '{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}'.format( \
         ab.data[0, 2, 67, 12], at.data[0, 3, 30, 25], abadv.data[0, -1, 47, 28],
         abpro.data[0, -1, 33, 17], brpror.data[0, 3, 28, 78], fladv.data[0, -1, 29, 29],
         fldif.data[0, 2, 12, 3], flpro.data[0, 1, 3, 4], jrpro3d.data[0, -1, 12, 10, 6],
         jrpro3d.data[1, 3, 10, 31, 12], jradvr.data[0, -1, 20, 11],
         jrdifr.data[0, 1, 15, 16], fldif.data_ic[0, -1, 81, 3], ab.data_ic[0, -1, 12, 12])

    # Write output for movie files
    file.write(st)
    file.close()


class TestdtBMovieOutputs(unittest.TestCase):

    def __init__(self, testName, dir, execCmd='mpirun -n 8 ../tmp/magic.exe'):
        super(TestdtBMovieOutputs, self).__init__(testName)
        self.dir = dir
        self.execCmd = execCmd
        self.startDir = os.getcwd()
        self.description = "Test dtB Movie outputs"

    def list2reason(self, exc_list):
        if exc_list and exc_list[-1][0] is self:
            return exc_list[-1][1]

    def setUp(self):
        # Cleaning when entering
        print('\nDirectory   :           {}'.format(self.dir))
        print('Description :           {}'.format(self.description))
        self.startTime = time.time()
        cleanDir(self.dir)
        for f in glob.glob('{}/*.start'.format(self.dir)):
            os.remove(f)
        os.chdir(self.dir)
        cmd = '{} {}/input.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

    def tearDown(self):
        # Cleaning when leaving
        os.chdir(self.startDir)
        cleanDir(self.dir)
        for f in glob.glob('{}/*.start'.format(self.dir)):
            os.remove(f)

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

    @unittest.skipUnless('MAGIC_HOME' in os.environ,
                         'MAGIC_HOME is not defined! source sourceme.sh!')
    def outputFileDiff(self):
        generateEkinFile('e_kin.test')
        datRef = readStack('{}/reference.out'.format(self.dir))
        datTmp = readStack('{}/e_kin.test'.format(self.dir))
        np.testing.assert_equal(datRef, datTmp)


if __name__ == '__main__':
    generateEkinFile('e_kin.test')
