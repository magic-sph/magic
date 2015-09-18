#
# If you use csh/tcsh:
#
# do a source sourceme.csh to initialize the correct PATHS for magic
#


if (! $?MAGIC_HOME) then
  unset _sourceme		# tabula rasa without MAGIC_HOME
  #
  # Try to identify position of the code's home directory:
  #
  foreach _dir ( . .. ../.. ../../.. ../../../.. magic)
    if ( (-e $_dir/sourceme.csh) && \
         (-d $_dir/src)          && \
	 (-d $_dir/bin)          && \
	 (-d $_dir/samples)         \
       ) then
      set back_again = `pwd`     
      cd $_dir; setenv MAGIC_HOME `pwd`; cd $back_again
      goto found
    endif
  end

  echo "sourceme.csh: Cannot locate home directory of MAGIC"
  echo "  Try sourcing me from the home directory itself, or set MAGIC_HOME"
  goto eof
endif
    
found:
if (! $?_sourceme_quiet) echo "MAGIC_HOME = <$MAGIC_HOME>"

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

./bin/autoConfigPython.sh

endif
