#
# If you use bash/zsh/ksh:
#
# do a source sourceme.sh to initialize the correct PATHS for magic
#

if [ -z $MAGIC_HOME ]; then
  unset _sourceme               # tabula rasa without MAGIC_HOME
  #
  # Try to identify position of the code's home directory:
  #
  for _dir in   . .. ../.. ../../.. ../../../.. \
	        MAGIC_mpi ; do
    if ( [ -e $_dir/sourceme.sh ] && \
         [ -d $_dir/src ]         && \
         [ -d $_dir/samples ]     && \
         [ -d $_dir/bin ]            \
       ); then
      unset cd   # some people are crazy enough to overload cd
      MAGIC_HOME=`cd $_dir; echo $PWD`
      export MAGIC_HOME
      break
    fi
  done
  unset _dir

  if [ -z $MAGIC_HOME ]; then # no success
    echo "sourceme.sh: Cannot locate home directory of MAGIC"
    echo " Try sourcing me from the home directory itself, or set MAGIC_HOME"
  fi
fi

if [ -z $_sourceme_quiet ]; then echo "MAGIC_HOME = <$MAGIC_HOME>"; fi

if [ -z $_sourceme ]; then	# called for the first time?
  # CDPATH="./:../:../../:../../../:$HOME"
  if [ -n $MAGIC_HOME ]; then

    #  Set shell path
    if [ -z $_sourceme_quiet ]; then echo "Adding $MAGIC_HOME to PATH"; fi
    PATH=${PATH}:$MAGIC_HOME/bin:$MAGIC_HOME/src

    #   Set PYTHONPATH
    if [ -z $PYTHONPATH ]; then
       PYTHONPATH="$MAGIC_HOME/python:$PWD/python"
    else 
       PYTHONPATH="$PYTHONPATH:$MAGIC_HOME/python:$PWD/python"
    fi
    # Remember that sourceme has been successfully run
    _sourceme="set"

    export PATH
    
  else
    if [ -n $MAGIC_HOME ]; then
      echo "Not adding $MAGIC_HOME/bin to PATH: not a directory"
    fi
  fi
fi
