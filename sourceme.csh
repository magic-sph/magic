#
# If you use csh/tcsh:
#
# do a source sourceme.csh to initialize the correct PATHS for magic
#
set _sourcecmd = ( $_ )

if ( ! $?MAGIC_HOME ) then
  
  unset _sourceme               # tabula rasa without MAGIC_HOME
  
  #
  # Try to identify position of the code's home directory:
  #
  
  set _dir = `readlink -f $_sourcecmd[2]`
  set _dir = `dirname $_dir` 
  
  if ( (-e $_dir/sourceme.csh) && \
       (-d $_dir/src)          && \
       (-d $_dir/bin)          && \
       (-d $_dir/samples)         \
     ) then
    setenv MAGIC_HOME "$_dir"
    echo $MAGIC_HOME
    unset _dir
  endif
endif

if (! $?_sourceme) then		# called for the fist time?
  if (-d $MAGIC_HOME/bin) then

    #  Set shell path
    if (! $?_sourceme_quiet) echo "Adding $MAGIC_HOME to PATH"
    set path = ( $path $MAGIC_HOME/bin $MAGIC_HOME/src)

    #  Set PYTHON path
    if ($?PYTHONPATH) then
      setenv PYTHONPATH "${PYTHONPATH}:${MAGIC_HOME}/python:${PWD}/python"
    else
      setenv PYTHONPATH "${MAGIC_HOME}/python:${PWD}/python"
    endif

    set _sourceme = 'set'

  else
    echo "Not adding $MAGIC_HOME/bin to PATH: not a directory"
  endif

  # For reading binary outputs
  setenv F_UFMTENDIAN "big"
  setenv GFORTRAN_CONVERT_UNIT "big_endian"

endif

if ( $?MAGIC_HOME ) then
  $MAGIC_HOME/bin/autoConfigPython.sh
endif
