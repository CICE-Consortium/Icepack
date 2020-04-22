#!/bin/sh

#-----------------------------------
# Generate documentation from icepack source code
# Will cut and paste all code between "autodocument_start" and "autodocument_end"
#   (defined by $startline and $endline below)
# Writes to source/user_guide/interfaces.include ($rstfile below)
# Expects $startline to include a section title string such as
#   ! autodocument_start icepack_subroutine_name
# That string (ie. icepack_subroutine_name) will be used as the section name 
#   in the rst documentation
# The pattern matching assumes the title is anything following the $startline string
#   and can include spaces
# The file interfaces.include is included automatically in the icepack documentation
# To update the documentation
#   cd to this directory
#   run this script (./generate_interfaces.sh)
#   git add the $rstfile
#   git commit, push, and create PR
# It is recommended that each section of icepack that is documented start with a
#   comment that includes a ! in the first chararacter postition.  This will insure
#   proper indenting of the section.
#
# For example:
#   ! autodocument_start test_section
#   ! this is an example
#
#      subroutine test(arg1)
#
#      integer, intent(in), optional :: arg1  ! test argument
#
#   ! autodocument_end
#
#-----------------------------------

startline="autodocument_start"
endline="autodocument_end"
inpfiles="../columnphysics/*.F90"
rstfile="./source/user_guide/interfaces.include"
dowrite=0

mv ${rstfile} ${rstfile}.orig

for file in ${inpfiles}; do

filename=`basename $file`
firstfileintfc=0

  while IFS= read -r line; do
    if [[ $line =~ .*$endline.* ]]; then
      dowrite=0
      echo " " >> $rstfile
    fi

    if [ $dowrite = 1 ]; then
      echo "  $line" >> $rstfile
    fi

    if [[ $line =~ .*$startline.* ]]; then
      dowrite=1
      title=`echo ${line} | sed "s/.*${startline}[[:space:]]*\([^[:space:]]*.*\)$/\1/"`
      if [[ $title =~ ^[[:space:]]*$ ]]; then
         mv ${rstfile}.orig ${rstfile}
         echo " "
         echo "ERROR in $file"
         echo "  $line"
         echo "requires a section title"
         echo " "
         exit -9
      fi
      echo "$filename $title"
      if [ $firstfileintfc = 0 ]; then
         firstfileintfc=1
         echo "" >> $rstfile
         echo "${filename}" >> $rstfile
         echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> $rstfile
      fi
      echo "" >> $rstfile
      echo ".. _${title}:" >> $rstfile
      echo "" >> $rstfile
      echo "${title}" >> $rstfile
      echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" >> $rstfile
      echo ".. code-block:: fortran" >> $rstfile
      echo "" >> $rstfile
    fi

  done < "$file"

done

