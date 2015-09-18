#!/bin/sh

#
# This script tries to fill the magic.cfg for you
#

# Figure out which f2py executable is available on your machine
whichf2py () {
  if hash f2py3 2>/dev/null; then
    local f2pyExec="f2py3";
  elif hash f2py2 2>/dev/null; then
    local f2pyExec="f2py2";
  elif hash f2py 2>/dev/null; then
    local f2pyExec="f2py";
  else
    local f2pyExec="NotFound";
  fi
  echo $f2pyExec
}

# Awk magics to find out what are the available fortran compilers
whichf2pycompiler () {
  local pattern=`$1 -c --help-fcompiler | awk '/Compilers available/{flag=0}flag;/compilers found/{flag=1}' | awk -F'[=| ]' '{print $4}'`

  if  [ -z "$pattern" ]; then
    selectedCompiler="NotFound";
    echo $selectedCompiler
  else
    arr=($pattern)
    for compiler in "${arr[@]}"; do
      if [ $compiler == "intelem" ]; then
        sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
        break
      elif [ $compiler == "intele" ]; then
        sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
        break
      elif [ $compiler == "intel" ]; then
        sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
        break
      elif [ $compiler == "gnu95" ]; then
        sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
      elif [ $compiler == "pg" ]; then
        sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
      elif [ $compiler == "nag" ]; then
        sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
      elif [ $compiler == "gnu" ]; then
        sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
      fi
    done
    echo $selectedCompiler
  fi

}

# Check whether you can build the fortran libraries
buildLibs () {

  f2pyExec=$(whichf2py)
  if [ $f2pyExec != "NotFound" ]; then
    sed -i "s/f2pyexec.*/f2pyexec = $f2pyExec/g" $MAGIC_HOME/python/magic/magic.cfg

    local selectedCompiler=$(whichf2pycompiler $f2pyExec)
    if [ $selectedCompiler == "NotFound" ]; then
      sed -i "s/buildLib.*/buildLib = False/g" $MAGIC_HOME/python/magic/magic.cfg
    else
      sed -i "s/buildLib.*/buildLib = True/g" $MAGIC_HOME/python/magic/magic.cfg
    fi
  fi
}

getBackend () {
  local backendValue=`python $MAGIC_HOME/bin/testBackend.py 2>&1`
  if  [ -n "$backendValue" ]; then
      sed -i "s/backend.*/backend = $backendValue/g" $MAGIC_HOME/python/magic/magic.cfg
  fi

}


# Check if magic.cfg exists: if yes don't do anything, else create it from default
# and fill it with the estimated values
if [ ! -f $MAGIC_HOME/python/magic/magic.cfg ]; then
  if [ -w $MAGIC_HOME/python/magic/magic.cfg.default ]; then
    cp $MAGIC_HOME/python/magic/magic.cfg.default $MAGIC_HOME/python/magic/magic.cfg
    buildLibs
    getBackend
  fi
fi
