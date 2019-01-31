#!/bin/bash

#
# This script tries to fill the magic.cfg for you
#

# Spinner (in case of long loop)
sp="/-\|"
sc=0
spin() {
   printf "\b${sp:sc++:1}"
   ((sc==${#sp})) && sc=0
}
endspin() {
   printf "\r%s\n" "$@"
}

# Determines where matplotlib is installed
hasPythonMplDarwin() {
  if hash python 2>/dev/null; then
    local cmd=`python $MAGIC_HOME/bin/testBackend.py 2> /dev/null`
    if [ -n "$cmd" ]; then
      local backendValue=$cmd;
    else
      local backendValue="NotFound";
    fi
  else
    local backendValue="NotFound";
  fi
  echo $backendValue

}

hasPython2Mpl() {
  if hash python2 2>/dev/null; then
    local cmd=`python2 $MAGIC_HOME/bin/testBackend.py 2> /dev/null`
    if [ -n "$cmd" ]; then
      local backendValue=$cmd;
    else
      local backendValue="NotFound";
    fi
  else
    local backendValue="NotFound";
  fi
  echo $backendValue
}

hasPython3Mpl() {
  if hash python3 2>/dev/null; then
    local cmd=`python3 $MAGIC_HOME/bin/testBackend.py 2> /dev/null`
    if [ -n "$cmd" ]; then
      local backendValue=$cmd;
    else
      local backendValue="NotFound";
    fi
  else
    local backendValue="NotFound";
  fi
  echo $backendValue
}

# Figure out which f2py executable is available on your machine
hasf2py2 () {
  if hash f2py2 2>/dev/null; then
    local f2pyExec="f2py2";
  elif hash f2py 2>/dev/null; then
    local f2pyExec="f2py";
  elif hash f2py2.7 2>/dev/null; then
    local f2pyExec="f2py2.7";
  else
    local f2pyExec="NotFound";
  fi
  echo $f2pyExec
}

whichPython () {
  local cmd=`python -c 'import sys; print(sys.version_info[0])'`
  local pythonVersion=$cmd;
  echo $pythonVersion
}

hasf2py3 () {
  if hash f2py3 2>/dev/null; then
    local f2pyExec="f2py3";
  else
    local f2pyExec="NotFound";
  fi
  echo $f2pyExec
}

# Set sed command depending on OS type

whichSed () {

  if [[ $OSTYPE == *"darwin"* ]]; then
      SED='sed -i -e'
  else
      SED='sed -i'
  fi
}

# Awk magics to find out what are the available fortran compilers
whichf2pycompiler () {
  local pattern=`$1 -c --help-fcompiler | awk '/Compilers available/{flag=0}flag;/compilers found/{flag=1}' | awk -F'[=| ]' '{print $4}'`

  if  [ -z "$pattern" ]; then
    local selectedCompiler="NotFound";
    echo $selectedCompiler
  else
    local arr=($pattern)
    for compiler in "${arr[@]}"; do
      if [ $compiler == "intelem" ]; then
        $SED "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
        break
      elif [ $compiler == "intele" ]; then
        $SED "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
        break
      elif [ $compiler == "intel" ]; then
        $SED "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
        break
      elif [ $compiler == "gnu95" ]; then
        $SED "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
      elif [ $compiler == "pg" ]; then
        $SED "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
      elif [ $compiler == "nag" ]; then
        $SED "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
      elif [ $compiler == "gnu" ]; then
        $SED "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
        local selectedCompiler=$compiler;
      fi
    done
    echo $selectedCompiler
  fi

}

# Check whether you can build the fortran libraries
buildLibs () {

  local pythonVersion=$(whichPython)
  local f2py2Exec=$(hasf2py2)
  local f2py3Exec=$(hasf2py3)

  if [ $f2py2Exec != "NotFound" ]; then
    if [ $f2py3Exec != "NotFound" ]; then
      if [ $pythonVersion  == "3" ]; then
	local f2pyExec=$f2py3Exec
      else
        local f2pyExec=$f2py2Exec
      fi
    else
      local f2pyExec=$f2py2Exec
    fi
  else
    if [ $f2py3Exec == "NotFound" ]; then
      echo "f2py was not found"
      echo "f2py can't be set in $MAGIC_HOME/python/magic/magic.cfg"
      local f2pyExec="NotFound"
    else
      local f2pyExec=$f2py3Exec
    fi
  fi

  if [ $f2pyExec != "NotFound" ]; then
    $SED "s/f2pyexec.*/f2pyexec = $f2pyExec/g" $MAGIC_HOME/python/magic/magic.cfg

    local selectedCompiler=$(whichf2pycompiler $f2pyExec)
    if [ $selectedCompiler == "NotFound" ]; then
      $SED "s/buildLib.*/buildLib = False/g" $MAGIC_HOME/python/magic/magic.cfg
    else
      $SED "s/buildLib.*/buildLib = True/g" $MAGIC_HOME/python/magic/magic.cfg
    fi
  fi
}

# Get matplotlib backend
getBackend () {
  
  local backend2=$(hasPython2Mpl);
  local backend3=$(hasPython3Mpl);
  local backendDarwin=$(hasPythonMplDarwin);


  if [[ $OSTYPE == *"darwin"* ]]; then
    if [ $backendDarwin == "NotFound" ]; then
      echo "matplotlib was not found"
      echo "the backend can't be set in $MAGIC_HOME/python/magic/magic.cfg"
      local backendValue="NotFound"
    else
      local backendValue=$backendDarwin
    fi
  else
    if [ $backend2 != "NotFound" ]; then
      if [ $backend3 != "NotFound" ]; then
        echo "matplotlib is installed for both python2 and python3"
        echo "python2 is selected"
      fi
      local backendValue=$backend2
    else 
      if [ $backend3 == "NotFound" ]; then
        echo "matplotlib was not found"
        echo "the backend can't be set in $MAGIC_HOME/python/magic/magic.cfg"
        local backendValue="NotFound"
      else
        local backendValue=$backend3
      fi
    fi
  fi 


  if [ $backendValue != "NotFound" ]; then
      $SED "s/backend.*/backend = $backendValue/g" $MAGIC_HOME/python/magic/magic.cfg
  fi
}

ifGWDG() {

  local host_name=$HOSTNAME
 
  if [[ $host_name == *"gwd"* ]]
  then
      $SED "s/ccompiler.*/ccompiler = intelem/g" $MAGIC_HOME/python/magic/magic.cfg
  fi

}

# Check if magic.cfg exists: if yes don't do anything, else create it from default
# and fill it with the estimated values
if [ ! -f $MAGIC_HOME/python/magic/magic.cfg ]; then
  if [ -w $MAGIC_HOME/python/magic/magic.cfg.default ]; then
    cp $MAGIC_HOME/python/magic/magic.cfg.default $MAGIC_HOME/python/magic/magic.cfg
    echo "Trying to setup your PATH in $MAGIC_HOME/python/magic/magic.cfg..."
    echo "It might take some seconds, but it's gonna happen only once."
    whichSed
    buildLibs
    getBackend
    ifGWDG
  fi
fi

if [[ -f $MAGIC_HOME/python/magic/magic.cfg-e ]];then
    rm $MAGIC_HOME/python/magic/magic.cfg-e
fi
