#!/bin/sh


# Figure out which f2py executable is available on your machine
whichf2py () {
	if hash f2py3 2>/dev/null; then
		local f2pyExec="f2py3";
	elif hash f2py2 2>/dev/null; then
		local f2pyExec="f2py2";
	elif hash f2py 2>/dev/null; then
 		local f2pyExec="f2py";
	else
 		local f2pyExec="Not Found";
	fi
	echo $f2pyExec
}

# Awk magics to find out what are the available fortran compilers
whichf2pycompiler () {
	local pattern=`$1 -c --help-fcompiler | awk '/Compilers available/{flag=0}flag;/compilers found/{flag=1}' | awk -F'[=| ]' '{print $4}'`

	local compiler=''

	if  [ -z "$pattern" ]; then
		compiler="Not Found";
		echo $compiler
	else
		arr=($pattern)
		for compiler in "${arr[@]}"; do
			if [ $compiler == "intelem" ]; then
				sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
				selectedCompiler=$compiler;
				break
			elif [ $compiler == "intele" ]; then
				sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
				selectedCompiler=$compiler;
				break
			elif [ $compiler == "intel" ]; then
				sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
				selectedCompiler=$compiler;
				break
			elif [ $compiler == "gnu95" ]; then
				sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
				selectedCompiler=$compiler;
			elif [ $compiler == "pg" ]; then
				sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
				selectedCompiler=$compiler;
			elif [ $compiler == "nag" ]; then
				sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
				selectedCompiler=$compiler;
			elif [ $compiler == "gnu" ]; then
				sed -i "s/fcompiler.*/fcompiler = $compiler/g" $MAGIC_HOME/python/magic/magic.cfg
				selectedCompiler=$compiler;
			fi
		done
		echo $selectedCompiler
	fi

}

# Check whether you can build the fortran libraries
buildLibs () {

	f2pyExec=$(whichf2py)
	if [ $f2pyExec != "Not Found" ]; then
		sed -i "s/f2pyexec.*/f2pyexec = $f2pyExec/g" $MAGIC_HOME/python/magic/magic.cfg

		selectedCompiler=$(whichf2pycompiler $f2pyExec)
		if [ $selectedCompiler == "Not Found" ]; then
			sed -i "s/buildLib.*/buildLib = False/g" $MAGIC_HOME/python/magic/magic.cfg

		else
			sed -i "s/buildLib.*/buildLib = True/g" $MAGIC_HOME/python/magic/magic.cfg

		fi
	fi
}



if [ -w $MAGIC_HOME/python/magic/magic.cfg.default ]; then
	cp $MAGIC_HOME/python/magic/magic.cfg.default $MAGIC_HOME/python/magic/magic.cfg
	buildLibs
fi
