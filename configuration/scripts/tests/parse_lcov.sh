#!/bin/sh

INFILE=total.info

file="ZZZ"
tfunc=0
thit=0
tfnd=0
dacnt=0
declare -a func
declare -a lines
declare -a linee
declare -a fhit
declare -a ffnd
echo " "

while read -r LINE
do
  case $LINE in
    SF:*)
      file="${LINE##*/}"
      ;;
    FN:*)
      if [[ "$LINE" =~ ^FN:(.*),(.*),.*_MOD_(.*)$ ]]; then
        func[$tfunc]="${BASH_REMATCH[3]}"
        lines[$tfunc]="${BASH_REMATCH[1]}"
        linee[$tfunc]="${BASH_REMATCH[2]}"
        fhit[$tfunc]=0
        ffnd[$tfunc]=0
        tfunc=$(($tfunc + 1))
      fi
      ;;
    DA:*)
      if [[ "$LINE" =~ ^DA:(.*),(.*)$ ]]; then
        lnum="${BASH_REMATCH[1]}"
        lhit="${BASH_REMATCH[2]}"
#        echo "tcx3a ${linee[$dacnt]} $lnum $lhit $dacnt"
        if [[ $lnum -lt ${lines[$dacnt]} || $lnum -gt ${linee[$dacnt]} ]]; then
          cnt=0
          while [ $cnt -lt $tfunc ]; do
            if [[ $lnum -ge ${lines[$cnt]} && $lnum -le ${linee[$cnt]} ]]; then
              dacnt=$cnt
            fi
            cnt=$(($cnt + 1))
          done
        fi
#        echo "tcx3b $lnum $lhit $dacnt"
        ffnd[$dacnt]=$((ffnd[$dacnt] + 1))
        if [[ $lhit -gt 0 ]]; then
          fhit[$dacnt]=$((fhit[$dacnt] + 1))
        fi
      fi
      ;;
    LF:*)
      fnd="${LINE##*:}"
      ;;
    LH:*)
      hit="${LINE##*:}"
      ;;
    end_of_record)
      thit=$(($thit + $hit))
      tfnd=$(($tfnd + $fnd))
      percent=`echo "scale=2;  $hit * 100. / $fnd" | bc`
      printf "%6.2f%% %-42s %5d %5d\n" $percent $file $hit $fnd
      cnt=0
      while [ $cnt -lt $tfunc ]; do
        percent=`echo "scale=2;  $((fhit[$cnt] * 100 / ffnd[$cnt]))" | bc`
#        printf "%6.2f%% %-40s %5d %5d\n" 0 ${func[$cnt]} ${lines[$cnt]} ${linee[$cnt]}
        printf "%12.2f%% %-36s %5d %5d\n" $percent ${func[$cnt]} ${fhit[$cnt]} ${ffnd[$cnt]}
        cnt=$(($cnt + 1))
      done
      file="ZZZ"
      tfunc=0
      dacnt=0
      ;;
    *)
      ;;
  esac

done < "$INFILE"

echo " "
#echo "tcx1 $thit $tfnd"
percent=`echo "scale=2;  $thit * 100. / $tfnd" | bc`
printf "%6.2f%% %-32s %10d %10d\n" $percent **TOTAL** $thit $tfnd
echo " "
