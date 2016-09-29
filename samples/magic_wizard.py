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
import dynamo_benchmark.unitTest
import varProps.unitTest
import doubleDiffusion.unitTest
import testRestart.unitTest
import testMapping.unitTest
import testTruncations.unitTest
import fluxPerturbation.unitTest
import isothermal_nrho3.unitTest
import dynamo_benchmark_condICrotIC.unitTest
import varCond.unitTest
import hydro_bench_anel.unitTest
import testCoeffOutputs.unitTest
import testRMSOutputs.unitTest
import testGraphMovieOutputs.unitTest
import testTOGeosOutputs.unitTest



def getParser():
    """
    Get script option
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--level', action='store', dest='test_level', type=int,
                        default=-1, help='Test level')
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
    parser.add_argument('--use-hdf5', action='store_true', dest='use_hdf5', 
                        default=False, help='Use HDF5 to store restart files')
    parser.add_argument('--use-precond', action='store', dest='use_precond', 
                        type=bool, default=True, 
                        help='Use matrix preconditioning')
    parser.add_argument('--nthreads', action='store', dest='nthreads', type=int,
                        default=1,
                        help='Specify the number of threads (hybrid version)')
    parser.add_argument('--nranks', action='store', dest='nranks', type=int,
                        default=4, help='Specify the number of MPI ranks')
    parser.add_argument('--mpicmd', action='store', dest='mpicmd', type=str,
                        default='mpirun', help='Specify the mpi executable')

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
    print("                         '-'   '-'      ")
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
        if args.use_hdf5:
            hdf5_opt = '-DUSE_HDF5=yes'
        else:
            hdf5_opt = '-DUSE_HDF5=no'
    else:
        mpi_opt = '-DUSE_MPI=no'
        omp_opt = '-DUSE_OMP=no'
        hdf5_opt = '-DUSE_HDF5=no'

    # Compilation
    cmd = 'cmake %s/.. %s %s %s %s %s %s %s' % (startdir, mpi_opt, build_type,
                                                precond_opt, mkl_opt, omp_opt, 
                                                shtns_opt, hdf5_opt)
    print('  '+cmd)
    print('\n')
    sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))


def compile():
    """
    Compile the code
    """
    cmd = 'make -j'
    print('  '+cmd)
    print('\n')
    sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
            stderr=open(os.devnull, 'wb'))


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
        print('  HDF5      : %r' % args.use_hdf5)
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
    shutil.rmtree(execDir)


def getSuite(startdir, cmd, precision, args):
    """
    Construct test suite
    """
    suite = unittest.TestSuite()

    if args.test_level in [-1, 0]:
        # Initial state of the Boussinesq benchmark (non-conducting IC)
        suite.addTest(dynamo_benchmark.unitTest.DynamoBenchmark('outputFileDiff',
                                                  '%s/dynamo_benchmark' \
                                                  % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Variable properties
        suite.addTest(varProps.unitTest.VariableProperties('outputFileDiff',
                                                  '%s/varProps' % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Saturated state of the Boussinesq benchmark (conducting IC)
        suite.addTest(boussBenchSat.unitTest.BoussinesqBenchmarkTest(
                                                  'outputFileDiff',
                                                  '%s/boussBenchSat' % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Double Diffusion
        suite.addTest(doubleDiffusion.unitTest.DoubleDiffusion('outputFileDiff',
                                                  '%s/doubleDiffusion'\
                                                  % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
    if args.test_level in [-1, 1]:
        # Test restart file capabilities
        suite.addTest(testRestart.unitTest.TestRestart('outputFileDiff',
                                                  '%s/testRestart' % startdir, 
                                                  execCmd=cmd,
                                                  use_hdf5=args.use_hdf5,
                                                  precision=precision))
        # Test truncations
        suite.addTest(testTruncations.unitTest.TestTruncations('outputFileDiff',
                                                  '%s/testTruncations'\
                                                  % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Test grid mapping capabilities
        suite.addTest(testMapping.unitTest.TestMapping('outputFileDiff',
                                                  '%s/testMapping' % startdir, 
                                                  execCmd=cmd,
                                                  use_hdf5=args.use_hdf5,
                                                  precision=precision))
        # Check standard time series outputs
        suite.addTest(testOutputs.unitTest.OutputTest('outputFileDiff',
                                                  '%s/testOutputs' % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Check radial profiles
        suite.addTest(testRadialOutputs.unitTest.RadialOutputTest(
                                                  'outputFileDiff',
                                                  '%s/testRadialOutputs'\
                                                  % startdir, 
                                                  execCmd=cmd,
                                                  precision=1e-7))
    if args.test_level in [-1, 2]:
        # Check the anelastic non-magnetic benchmark
        suite.addTest(hydro_bench_anel.unitTest.AnelasticBenchmark(
                                                  'outputFileDiff',
                                                  '%s/hydro_bench_anel'\
                                                  % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
    if args.test_level in [-1, 3]:
        # Check heat flux pattern in a Boussinesq model
        suite.addTest(fluxPerturbation.unitTest.HeatFluxPattern(
                                                  'outputFileDiff',
                                                  '%s/fluxPerturbation'\
                                                  % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Check an anelastic case with a zero Gruneisen (isothermal) parameter
        suite.addTest(isothermal_nrho3.unitTest.ZeroGruneisen(
                                                  'outputFileDiff',
                                                  '%s/isothermal_nrho3'\
                                                  % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Check conducting and rotating IC
        suite.addTest(dynamo_benchmark_condICrotIC.unitTest.ConductingRotatingIC(
                                                  'outputFileDiff',
                                                  '%s/dynamo_benchmark_condICrotIC'\
                                                  % startdir,
                                                  execCmd=cmd,
                                                  precision=precision))
        # Check variable electrical conductivity
        suite.addTest(varCond.unitTest.VariableConductivity('outputFileDiff',
                                                  '%s/varCond' % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
    if args.test_level in [-1, 4]:
        # Check CMB and r coefficients
        suite.addTest(testCoeffOutputs.unitTest.TestCoeffOutputs(
                                                  'outputFileDiff',
                                                  '%s/testCoeffOutputs'\
                                                  % startdir, 
                                                  execCmd=cmd))
        # Check RMS force balance
        suite.addTest(testRMSOutputs.unitTest.TestRMSOutputs('outputFileDiff',
                                                  '%s/testRMSOutputs'\
                                                  % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Check Graphic and Movie outputs
        suite.addTest(testGraphMovieOutputs.unitTest.TestGraphicMovieOutputs(
                                                  'outputFileDiff',
                                                  '%s/testGraphMovieOutputs' \
                                                  % startdir, 
                                                  execCmd=cmd))
        # Check TO and Geos outputs
        suite.addTest(testTOGeosOutputs.unitTest.TestTOGeosOutputs(
                                                  'outputFileDiff',
                                                  '%s/testTOGeosOutputs' \
                                                  % startdir, 
                                                  execCmd=cmd))

    return suite


if __name__ == '__main__':
    precision = 1e-8 # relative tolerance between expected and actual result
    startdir = os.getcwd()
    execDir = '%s/tmp' % startdir # where MagIC will be built

    parser = getParser()
    args = parser.parse_args()

    # Initialisation
    wizard()

    # Display environment variables
    print('0.    Environment     ')
    print('----------------------')
    get_env(args)

    # Run cmake
    print('1. cmake configuration')
    print('----------------------')
    cmake(args, startdir, execDir)

    # Compile the code
    print('2.    compilation     ')
    print('----------------------')
    compile()

    # Determine the execution command
    cmd = get_exec_cmd(args, execDir)

    # Run the auto-test suite
    print('3.   Auto-tests       ')
    print('----------------------')
    suite = getSuite(startdir, cmd, precision, args)
    runner = unittest.TextTestRunner(verbosity=0)
    runner.run(suite)

    # Clean build directory
    clean_exec_dir(execDir)
