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
    for f in glob.glob('{}/*.start'.format(dir)):
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
    from magic import MagicGraph, Movie

    gr = MagicGraph(ivar=1)

    file = open(fileName, 'w')
    st = '{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}'.format( \
         gr.entropy[10, 13, 3], gr.Br[33, 0, 2], gr.Btheta[3, 11, 11],
         gr.Bphi[34, 12, 25], gr.vr[11, 15, 2], gr.vtheta[33, 33, 3],
         gr.vphi[32, 33, 7], gr.Br_ic[5, 101, 3], gr.Btheta_ic[40, 39, 13] )

    #Write output for graphic files
    file.write(st+'\n')

    av = Movie(file='AV_mov.start', iplot=False)
    ahf = Movie(file='AHF_mov.start', iplot=False)
    brcmb = Movie(file='Br_CMB_mov.start', iplot=False)
    vtr = Movie(file='Vt_R=1_338461_mov.start', iplot=False)
    vreq = Movie(file='Vr_EQU_mov.start', iplot=False)
    teq = Movie(file='T_EQU_mov.start', iplot=False)
    vortz = Movie(file='VorZ_EQU_mov.start', iplot=False)
    hel = Movie(file='HE_R=1_438461_mov.start', iplot=False)
    Bteq = Movie(file='Bt_EQU_mov.start', iplot=False)
    Br3D = Movie(file='Br_3D_mov.start', iplot=False)
    vt3D = Movie(file='Vt_3D_mov.start', iplot=False)
    fl = Movie(file='FL_mov.start', iplot=False)
    vr_slice = Movie(file='Vr_P=0_000000_mov.start', iplot=False)
    vp_slice = Movie(file='Vp_P=0_000000_mov.start', iplot=False)
    Br_slice = Movie(file='Br_P=0_000000_mov.start', iplot=False)

    st = '{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}'.format( \
         av.data[0, 1, 121, 12], ahf.data[0, 0, 99, 33], brcmb.data[0, 1, 47, 128],
         vtr.data[0, 0, 33, 87], vreq.data[0, 0, 28, 31], teq.data[0, 1, 29, 29],
         vortz.data[0, 1, 1, 1], hel.data[0, 0, 3, 4], Bteq.data_ic[0, -1, 10, 6],
         Br3D.data[0, -1, 10, 31, 12], Br3D.data_ic[0, -1, 20, 11, 5],
         vt3D.data[0, 0, 15, 16, 13], fl.data[0, -1, 10, 12],
         vr_slice.data[0, -1, 10, 10], vp_slice.data[0, -1, 73, 13],
         Br_slice.data[0, 0, 212, 13], Br_slice.data_ic[0, -1, 193, 8])

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
        print('\nDirectory   :           {}'.format(self.dir))
        print('Description :           {}'.format(self.description))
        self.startTime = time.time()
        cleanDir(self.dir)
        os.chdir(self.dir)
        cmd = '{} {}/input.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

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

    @unittest.skipUnless('MAGIC_HOME' in os.environ,
                         'MAGIC_HOME is not defined! source sourceme.sh!')
    def outputFileDiff(self):
        generateEkinFile('e_kin.test')
        datRef = readStack('{}/reference.out'.format(self.dir))
        datTmp = readStack('{}/e_kin.test'.format(self.dir))
        np.testing.assert_equal(datRef, datTmp)


if __name__ == '__main__':
    generateEkinFile('e_kin.test')
