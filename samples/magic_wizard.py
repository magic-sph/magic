#!/usr/bin/env python
from __future__ import print_function
import argparse
import os
import shutil
import sys
import subprocess as sp
import numpy as np
import unittest
import testOutputs.unitTest
import testRadialOutputs.unitTest
import boussBenchSat.unitTest
import precession.unitTest
import dynamo_benchmark.unitTest
import full_sphere.unitTest
import varProps.unitTest
import finite_differences.unitTest
import time_schemes.unitTest
import doubleDiffusion.unitTest
import testRestart.unitTest
import testMapping.unitTest
import testTruncations.unitTest
import fluxPerturbation.unitTest
import isothermal_nrho3.unitTest
import dynamo_benchmark_condICrotIC.unitTest
import varCond.unitTest
import hydro_bench_anel.unitTest
import couetteAxi.unitTest
import testCoeffOutputs.unitTest
import testRMSOutputs.unitTest
import testGraphMovieOutputs.unitTest
import testTOGeosOutputs.unitTest

__version__ = '1.0.1'


def getParser():
    """
    Get script option
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s '+__version__, 
                        help="Show program's version number and exit.")
    parser.add_argument('--level', action='store', dest='test_level', type=int,
                        default=-1, help='Test level, use -2 for more info')
    parser.add_argument('--use-debug-flags', action='store_true', 
                        dest='use_debug_flags', 
                        default=False, help='Use compilation debug flags')
    parser.add_argument('--use-mpi', action='store_true', dest='use_mpi', 
                        default=False, help='Use MPI')
    parser.add_argument('--use-openmp', action='store_true', dest='use_openmp', 
                        default=False, help='Use the hybrid version')
    parser.add_argument('--use-mkl', action='store_true', dest='use_mkl',   
                        default=False, 
                        help='Use the MKL for FFTs and Lapack calls')
    parser.add_argument('--use-shtns', action='store_true', dest='use_shtns', 
                        default=False, help='Use SHTns for Legendre transforms')
    parser.add_argument('--use-precond', action='store', dest='use_precond', 
                        type=bool, default=True, 
                        help='Use matrix preconditioning')
    parser.add_argument('--nranks', action='store', dest='nranks', type=int,
                        default=4, help='Specify the number of MPI ranks')
    parser.add_argument('--nthreads', action='store', dest='nthreads', type=int,
                        default=1,
                        help='Specify the number of threads (hybrid version)')
    parser.add_argument('--mpicmd', action='store', dest='mpicmd', type=str,
                        default='mpirun', help='Specify the mpi executable')
    parser.add_argument('--name', action='store', dest='test_name', type=str,
                        default='all', help='Runs only the specified tests')
    parser.add_argument('--use-exec', action='store', dest='use_exec', type=str,
                        default='', help='Uses executable from <path> instead of compiling it from source')
    parser.add_argument('--extra-arg', action='store', dest='extra_arg', type=str,
                        default='', 
                        help='extra command line arguments to be passed to magic.exe')
    parser.add_argument('--log', action='store_true', dest='log', 
                        default=False, 
                        help='keeps the output for all runs')

    return parser


def wizard():
    print('\n')
    print("MagIC's auto-tests suite")
    print('\n')
    print("                      /\                ")
    print("                     /  \               ")
    print("                    |    |              ")
    print("                  --:'''':--            ")
    print("                    :'_' :              ")
    print('                    _:"":\___           ')
    print("     ' '      ____.' :::     '._        ")
    print("    . *=====<<=)           \    :       ")
    print("     .  '      '-'-'\_      /'._.'      ")
    print('                      \====:_ ""        ')
    print("                     .'     \\          ")
    print("                    :       :           ")
    print("                   /   :    \           ")
    print("                  :   .      '.         ")
    print("                  :  : :      :         ")
    print("                  :__:-:__.;--'         ")
    print("                   '-'   '-'            ")
    print('\n')


def cmake(args, startdir, execDir):
    """
    Run cmake
    """
    if not os.path.exists(execDir):
        os.mkdir(execDir)
    else:
        shutil.rmtree(execDir)
        os.mkdir(execDir)
    os.chdir(execDir)

    if args.use_debug_flags:
         build_type='-DCMAKE_BUILD_TYPE=Debug'
    else:
         build_type='-DCMAKE_BUILD_TYPE=Release'

    if args.use_precond:
        precond_opt = '-DUSE_PRECOND=yes'
    else:
        precond_opt = '-DUSE_PRECOND=no'

    if args.use_mkl:
        mkl_opt = '-DUSE_FFTLIB=MKL -DUSE_LAPACKLIB=MKL'
    else:
        mkl_opt = '-DUSE_FFTLIB=JW -DUSE_LAPACKLIB=JW'

    if args.use_shtns:
        shtns_opt = '-DUSE_SHTNS=yes'
    else:
        shtns_opt = '-DUSE_SHTNS=no'

    if args.use_mpi:
        mpi_opt = '-DUSE_MPI=yes'
        if args.use_openmp:
            omp_opt = '-DUSE_OMP=yes'
        else:
            omp_opt = '-DUSE_OMP=no'
    else:
        mpi_opt = '-DUSE_MPI=no'
        omp_opt = '-DUSE_OMP=no'

    # Compilation
    cmd = 'cmake %s/.. %s %s %s %s %s %s' % (startdir, mpi_opt, build_type,
                                                precond_opt, mkl_opt, omp_opt, 
                                                shtns_opt)
    print('  '+cmd)
    print('\n')
    sp.call(cmd, shell=True, stdout=args.cmake_outFile, stderr=sp.STDOUT)
    args.cmake_outFile.close()


def compile(args):
    """
    Compile the code
    """
    cmd = 'make -j'
    if args.log: cmd += ' VERBOSE=1'
    print('  '+cmd)
    print('\n')
    sp.call(cmd, shell=True, stdout=args.make_outFile, stderr=sp.STDOUT)
    args.make_outFile.close()


def get_env(args):
    """
    Display the environment variables
    """
    if 'FC' in os.environ:
        fortran_comp = os.environ['FC']
    else:
        fortran_comp = 'FC is not defined. Default Fortran compiler will be used!'
    if 'CC' in os.environ:
        c_comp = os.environ['CC']
    else:
        c_comp = 'CC is not defined. Default C compiler will be used!'
        
    print('  FC        : %s' % fortran_comp)
    print('  CC        : %s' % c_comp)
    print('  MPI       : %r' % args.use_mpi)
    if args.use_mpi:
        print('  nranks    : %i' % args.nranks)
        print('  mpi exec  : %s' % args.mpicmd)
        print('  OPENMP    : %r' % args.use_openmp)
        if args.use_openmp:
            print('  nThreads  : %i' % args.nthreads)
        print('  MKL       : %r' % args.use_mkl)
        print('  SHTNS     : %r' % args.use_shtns)
    print('\n')


def get_exec_cmd(args, execDir):
    """
    Determine execution command
    """
    magicExec = '%s/magic.exe' % execDir

    if args.use_mpi: # With MPI
        if args.use_openmp:
            os.environ['OMP_NUM_THREADS'] = str(args.nthreads)
            os.environ['KMP_STACKSIZE'] = '1g'
            os.environ['KMP_AFFINITY'] = 'noverbose,granularity=core,compact'
            os.environ['MP_BINDPROC'] = 'no'
            os.environ['I_MPI_PIN_PROCESSOR_LIST'] = 'socket'
        else:
            os.environ['I_MPI_PIN_PROCESSOR_LIST'] = 'allcores'

        execCmd = '%s -n %i %s' % (args.mpicmd, args.nranks, magicExec)
    else: # Without MPI
        execCmd = '%s' % magicExec

    return execCmd


def clean_exec_dir(execDir):
    """
    Remove build directory
    """
    #shutil.rmtree(execDir)


def getSuite(startdir, cmd, precision, args):
    """
    Construct test suite
    """
    suite = unittest.TestSuite()
    
    test_name = args.test_name
    if args.test_name == 'all':
       test_level = args.test_level
    else:
       test_level = -1

    if test_level in [-1, 0]:
        # Initial state of the Boussinesq benchmark (non-conducting IC)
        if test_name == 'all' or 'dynamo_benchmark' in test_name:
           suite.addTest(dynamo_benchmark.unitTest.DynamoBenchmark('outputFileDiff',
                                                     '%s/dynamo_benchmark' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
        # Variable properties
        if test_name == 'all' or 'varProps' in test_name:
           suite.addTest(varProps.unitTest.VariableProperties('outputFileDiff',
                                                     '%s/varProps' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
           
        # Finite differences
        if test_name == 'all' or 'finite_differences' in test_name:
           suite.addTest(finite_differences.unitTest.FiniteDifferences('outputFileDiff',
                                                     '%s/finite_differences' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
           

        # Saturated state of the Boussinesq benchmark (conducting IC)
        if test_name == 'all' or 'boussBenchSat' in test_name:
           suite.addTest(boussBenchSat.unitTest.BoussinesqBenchmarkTest(
                                                     'outputFileDiff',
                                                     '%s/boussBenchSat' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
        # Double Diffusion
        if test_name == 'all' or 'doubleDiffusion' in test_name:
           suite.addTest(doubleDiffusion.unitTest.DoubleDiffusion('outputFileDiff',
                                                     '%s/doubleDiffusion' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
        # Axisymmetric run (spherical Couette)
        if test_name == 'all' or 'couetteAxi' in test_name:
           suite.addTest(couetteAxi.unitTest.CouetteAxi('outputFileDiff',
                                                     '%s/couetteAxi'\
                                                     % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
        # Precession with Po=-0.01
        if test_name == 'all' or 'outputFileDiff' in test_name:
          suite.addTest(precession.unitTest.PrecessionTest(
                                                    'outputFileDiff',
                                                    '%s/precession' % startdir, 
                                                    execCmd=cmd, log=args.log, 
                                                    precision=precision))

        # First full sphere benchmark from Marti et al. (2014)
        if test_name == 'all' or 'full_sphere' in test_name:
          if args.use_shtns:
              suite.addTest(full_sphere.unitTest.FullSphere('outputFileDiff',
                                                        '%s/full_sphere' \
                                                        % startdir, 
                                                        execCmd=cmd, log=args.log, 
                                                        precision=precision))
    if args.test_level in [-1, 1]:
        # Test restart file capabilities
        if test_name == 'all' or 'testRestart' in test_name:
           suite.addTest(testRestart.unitTest.TestRestart('outputFileDiff',
                                                     '%s/testRestart' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
        # Test truncations
        if test_name == 'all' or 'testTruncations' in test_name:
           suite.addTest(testTruncations.unitTest.TestTruncations('outputFileDiff',
                                                     '%s/testTruncations' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
        # Test grid mapping capabilities
        if test_name == 'all' or 'testMapping' in test_name:
           suite.addTest(testMapping.unitTest.TestMapping('outputFileDiff',
                                                     '%s/testMapping' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
        # Check standard time series outputs
        if test_name == 'all' or 'testOutputs' in test_name:
           suite.addTest(testOutputs.unitTest.OutputTest('outputFileDiff',
                                                     '%s/testOutputs' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
        # Check radial profiles
        if test_name == 'all' or 'testRadialOutputs' in test_name:
           suite.addTest(testRadialOutputs.unitTest.RadialOutputTest(
                                                     'outputFileDiff',
                                                     '%s/testRadialOutputs' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=1e-7))
    if test_level in [-1, 2]:
        # Check the anelastic non-magnetic benchmark
        if test_name == 'all' or 'hydro_bench_anel' in test_name:
           suite.addTest(hydro_bench_anel.unitTest.AnelasticBenchmark(
                                                     'outputFileDiff',
                                                     '%s/hydro_bench_anel' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
    if test_level in [-1, 3]:
        # Check heat flux pattern in a Boussinesq model
        if test_name == 'all' or 'fluxPerturbation' in test_name:
           suite.addTest(fluxPerturbation.unitTest.HeatFluxPattern(
                                                     'outputFileDiff',
                                                     '%s/fluxPerturbation' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
        # Check an anelastic case with a zero Gruneisen (isothermal) parameter
        if test_name == 'all' or 'isothermal_nrho3' in test_name:
           suite.addTest(isothermal_nrho3.unitTest.ZeroGruneisen(
                                                     'outputFileDiff',
                                                     '%s/isothermal_nrho3' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
        # Check conducting and rotating IC
        if test_name == 'all' or 'dynamo_benchmark_condICrotIC' in test_name:
           suite.addTest(dynamo_benchmark_condICrotIC.unitTest.ConductingRotatingIC(
                                                     'outputFileDiff',
                                                     '%s/dynamo_benchmark_condICrotIC' % startdir,
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
        # Check variable electrical conductivity
        if test_name == 'all' or 'varCond' in test_name:
           suite.addTest(varCond.unitTest.VariableConductivity('outputFileDiff',
                                                     '%s/varCond' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
    if test_level in [-1, 4]:
        # Check CMB and r coefficients
        if test_name == 'all' or 'testCoeffOutputs' in test_name:
           suite.addTest(testCoeffOutputs.unitTest.TestCoeffOutputs(
                                                     'outputFileDiff',
                                                     '%s/testCoeffOutputs' % startdir, 
                                                     execCmd=cmd))
        # Check RMS force balance
        if test_name == 'all' or 'testRMSOutputs' in test_name:
           suite.addTest(testRMSOutputs.unitTest.TestRMSOutputs('outputFileDiff',
                                                     '%s/testRMSOutputs' % startdir, 
                                                     execCmd=cmd, log=args.log, 
                                                     precision=precision))
        # Check Graphic and Movie outputs
        if test_name == 'all' or 'testGraphMovieOutputs' in test_name:
           suite.addTest(testGraphMovieOutputs.unitTest.TestGraphicMovieOutputs(
                                                     'outputFileDiff',
                                                     '%s/testGraphMovieOutputs' % startdir, 
                                                     execCmd=cmd))
        # Check TO and Geos outputs
        if test_name == 'all' or 'testTOGeosOutputs' in test_name:
           suite.addTest(testTOGeosOutputs.unitTest.TestTOGeosOutputs(
                                                     'outputFileDiff',
                                                     '%s/testTOGeosOutputs' % startdir, 
                                                     execCmd=cmd))

    return suite

def printLevelInfo():

        print("Level                               Test Description                               ")
        print("                                                                                   ")

        print("  0                 Boussinesq benchmark (start), non-conducting IC                ")
        print("                    Variable transport properties (anelastic, both Cheb and FD)    ")
        print("                    Test finite differences (restart from Cheb)                    ")
        print("                    Boussinesq benchmark: saturated state                          ")
        print("                    Double-diffusive convection (Breuer et al.)                    ")
        print("                    Test running an axisymmetric Couette flow                      ")
        print("                    Test precession with Po=-0.01 and E=1e-3                       ")
        print("                    Test full sphere benchmark from Marti et al. (2014)            ")
        print("                                                                                   ")
        print("                                                                                   ")


        print("  1                 Test restarting from a check point                             ")
        print("                    Test various truncations                                       ")
        print("                    Test the remapping of the grid at restart                      ")
        print("                    Test classical time-series outputs                             ")
        print("                    Test time-averaged radial profiles                             ")
        print("                                                                                   ")
        print("                                                                                   ")


        print("  2                 Anelastic non-magnetic benchmark (Jones et al.)                ")
        print("                                                                                   ")
        print("                                                                                   ")


        print("  3                 Heat flux pattern at the outer boundary                        ")
        print("                    Anelastic with zero Gruneisen (isothermal)                     ")
        print("                    Conducting and rotating inner core                             ")
        print("                    Anelastic model with variable conductivity                     ")
        print("                                                                                   ")
        print("                                                                                   ")


        print("  4                 Test coeff outputs                                             ")
        print("                    Test r.m.s. force balance calculation                          ")
        print("                    Test Graphic and Movie outputs                                 ")
        print("                    Test TO and Geos outputs                                       ")


def define_outputs(args):
    """
    Defines output for cmake and make
    """
    if args.log:
        args.cmake_outFile = open('./cmake.out','w')
        args.make_outFile = open('./make.out','w')
    else:
        args.cmake_outFile = open(os.devnull,'wb')
        args.make_outFile = open(os.devnull,'wb')
        
    if args.use_exec == '':
        args.exec_path = ''
    else:
        args.exec_path = os.path.abspath(startdir + '/' + args.use_exec + '/magic.exe')
    
    return args


if __name__ == '__main__':
    precision = 1e-8 # relative tolerance between expected and actual result
    startdir = os.getcwd()
    execDir = '%s/tmp' % startdir # where MagIC will be built

    parser = getParser()
    args = parser.parse_args()
    args = define_outputs(args)

    # Initialisation
    wizard()

    # Display environment variables
    print('0.    Environment     ')
    print('----------------------')
    get_env(args)

    if args.test_level == -2:
        printLevelInfo()
        sys.exit()

    if args.use_exec != '':
      # Symbolic linking instead of compiling
      print('1/2. Symbolic Linking   ')
      print('----------------------')
      print('Source: ' + os.path.abspath(args.exec_path))
      print('Destination: ' + os.path.abspath(execDir+'/magic.exe\n'))
      try:
         os.mkdir(execDir)
      except:
         pass
      try:
         os.unlink(execDir+'/magic.exe')
      except:
         pass
      os.symlink( args.exec_path , execDir+'/magic.exe' )
    else:
      # Run cmake
      print('1. cmake configuration')
      print('----------------------')
      cmake(args, startdir, execDir)

      # Compile the code
      print('2.    compilation     ')
      print('----------------------')
      compile(args)

    # Determine the execution command
    cmd = get_exec_cmd(args, execDir)

    # Run the auto-test suite
    print('3.   Auto-tests       ')
    print('----------------------')
    suite = getSuite(startdir, cmd, precision, args)
    runner = unittest.TextTestRunner(verbosity=0)
    ret = not runner.run(suite).wasSuccessful()

    # Clean build directory
    clean_exec_dir(execDir)

    sys.exit(ret)
