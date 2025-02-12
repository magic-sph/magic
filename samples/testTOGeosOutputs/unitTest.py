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
        dat = np.asarray(cut, np.float64)
        #out.append(dat)
        out = np.append(out, dat)
    return out


def generateEkinFile(fileName='e_kin.test'):
    from magic import TOMovie, MagicTOHemi, Movie

    # Write output for TO_mov file
    to = TOMovie(file='TO_mov.start', iplot=False)
    out = 'tmp'
    file = open(out, 'w')
    st = '{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}'.format( to.asVphi[0, 13, 3],
         to.rey[1, 21, 22], to.adv[1, 52, 11], to.visc[0, 12, 25],
         to.lorentz[0, 73, 30], to.coriolis[1, 33, 3], to.dtVp[1, 88, 7] )
    file.write(st+'\n')

    # TOnhs.TAG
    to = MagicTOHemi(hemi='n')
    st = '{:.4f} {:.4f} {:.4f} {:.4f}'.format(to.vp[2, 18], to.rstr[3, 11],
                                              to.astr[0, 30], to.LF[2, 21])
    file.write(st+'\n')

    # TOshs.TAG
    to = MagicTOHemi(hemi='s')
    st = '{:.4f} {:.4f} {:.4f} {:.4f}'.format(to.dvp[3, 12], to.viscstr[1, 9],
                                              to.tay[4, 27], to.vpr[2, 21])
    file.write(st+'\n')

    # Geos movies
    vs_geos = Movie(file='geosVS_mov.start', iplot=False)
    vp_geos = Movie(file='geosVPHI_mov.start', iplot=False)
    vortz_geos = Movie(file='geosVorZ_mov.start', iplot=False)

    st = '{:.4f} {:.4f} {:.4f}'.format(vs_geos.data[0, -1, 24, 12],
                                       vp_geos.data[0, -1, 37, 32],
                                       vortz_geos.data[0, -1, 7, 5])
    file.write(st+'\n')
    file.close()

    # Cat e_kin.test + misc
    with open(fileName, 'w') as outFile:
        sp.call(['cat', 'geos.start', 'Tay.start', out], stdout=outFile)

    os.remove('tmp')


class TestTOGeosOutputs(unittest.TestCase):

    def __init__(self, testName, dir, execCmd='mpirun -n 8 ../tmp/magic.exe'):
        super(TestTOGeosOutputs, self).__init__(testName)
        self.dir = dir
        self.execCmd = execCmd
        self.startDir = os.getcwd()
        self.description = "Test TO and Geos outputs"

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
