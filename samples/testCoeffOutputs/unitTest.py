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
    from magic import MagicCoeffCmb, MagicCoeffR, MagicPotential

    file = open(fileName, 'w')

    # Classical CMB file
    cmb = MagicCoeffCmb(tag='start', iplot=False, quiet=True)
    st = '{:.5e} {:.5e} {:.5e} {:.5e} {:.5e}'.format( cmb.glm[0, cmb.idx[1, 0]],
                                                      cmb.glm[2, cmb.idx[1, 0]],
                                                      cmb.hlm[-1, cmb.idx[5, 0]],
                                                      cmb.glm[2, cmb.idx[5, 4]],
                                                      cmb.hlm[3, cmb.idx[5, 4]] )
    file.write(st+'\n')

    # Secular variation CMB file
    cmb = MagicCoeffCmb(tag='start', iplot=False, sv=True, quiet=True)
    st = '{:.5e} {:.5e} {:.5e} {:.5e} {:.5e}'.format( cmb.glm[1, cmb.idx[1, 0]],
                                                      cmb.glm[2, cmb.idx[1, 0]],
                                                      cmb.hlm[-1, cmb.idx[5, 0]],
                                                      cmb.glm[2, cmb.idx[5, 4]],
                                                      cmb.hlm[3, cmb.idx[5, 4]] )
    file.write(st+'\n')

    # Time-averaged CMB file
    cmb = MagicCoeffCmb(tag='start',  iplot=False, ave=True, quiet=True)
    st = '{:.5e} {:.5e} {:.5e} {:.5e}'.format( cmb.glm[0, cmb.idx[1, 0]],
                                               cmb.hlm[0, cmb.idx[5, 0]],
                                               cmb.glm[0, cmb.idx[5, 4]],
                                               cmb.hlm[0, cmb.idx[5, 4]] )
    file.write(st+'\n')

    # Coeff at depth
    coeff = MagicCoeffR(tag='start',  field='V', r=1, iplot=False, quiet=True)
    st = '{:.5e} {:.5e} {:.5e} {:.5e}'.format( coeff.wlm[0, coeff.idx[4, 4]].real,
                                               coeff.wlm[0, coeff.idx[4, 4]].imag,
                                               coeff.zlm[0, coeff.idx[1, 0]].real,
                                               coeff.wlm[0, coeff.idx[2, 0]].imag )
    file.write(st+'\n')

    coeff = MagicCoeffR(tag='start',  field='V', r=2, iplot=False, quiet=True)
    st = '{:.5e} {:.5e} {:.5e} {:.5e}'.format( coeff.wlm[0, coeff.idx[4, 4]].real,
                                               coeff.wlm[0, coeff.idx[4, 4]].imag,
                                               coeff.zlm[0, coeff.idx[1, 0]].real,
                                               coeff.wlm[0, coeff.idx[2, 0]].imag )
    file.write(st+'\n')

    coeff = MagicCoeffR(tag='start',  field='B', r=1, iplot=False, quiet=True)
    st = '{:.5e} {:.5e} {:.5e} {:.5e}'.format( coeff.wlm[0, coeff.idx[1, 0]].real,
                                               coeff.wlm[0, coeff.idx[1, 0]].imag,
                                               coeff.zlm[0, coeff.idx[5, 0]].real,
                                               coeff.wlm[0, coeff.idx[5, 4]].imag )
    file.write(st+'\n')

    coeff = MagicCoeffR(tag='start',  field='B', r=2, iplot=False, quiet=True)
    st = '{:.5e} {:.5e} {:.5e} {:.5e}'.format( coeff.wlm[0, coeff.idx[1, 0]].real,
                                               coeff.wlm[0, coeff.idx[1, 0]].imag,
                                               coeff.zlm[0, coeff.idx[4, 0]].real,
                                               coeff.wlm[0, coeff.idx[5, 4]].imag )
    file.write(st+'\n')

    pot = MagicPotential(field='B', verbose=False)
    st = '{:.5e} {:.5e} {:.5e} {:.5e}'.format( pot.pol[5, 3].real,
                                               pot.pol[127, 15].imag,
                                               pot.tor[183, 16].real,
                                               pot.tor[121, 9].imag )
    file.write(st+'\n')

    pot = MagicPotential(field='V', verbose=False)
    st = '{:.5e} {:.5e} {:.5e} {:.5e}'.format( pot.pol[65, 15].real,
                                               pot.pol[125, 17].imag,
                                               pot.tor[63, 12].real,
                                               pot.tor[184, 13].imag )
    file.write(st+'\n')

    pot = MagicPotential(field='T', verbose=False)
    st = '{:.5e} {:.5e}'.format( pot.pol[0, 15].real, pot.pol[65, 9].real )
    file.write(st+'\n')

    file.close()


class TestCoeffOutputs(unittest.TestCase):

    def __init__(self, testName, dir, execCmd='mpirun -n 8 ../tmp/magic.exe'):
        super(TestCoeffOutputs, self).__init__(testName)
        self.dir = dir
        self.execCmd = execCmd
        self.startDir = os.getcwd()
        self.description = "Test coeff outputs"

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
