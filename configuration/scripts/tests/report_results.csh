#!/bin/csh -f

set wikirepo = "https://github.com/CICE-Consortium/Test-Results.wiki.git"
set wikiname = Test-Results.wiki

set tsubdir = icepack_testing
set hfile = "icepack_by_hash"
set mfile = "icepack_by_mach"
set vfile = "icepack_by_vers"
set bfile = "icepack_by_bran"

rm -r -f ${wikiname}
git clone ${wikirepo} ${wikiname}

set repo = `grep "#repo = " results.log | cut -c 9-`
set bran = `grep "#bran = " results.log | cut -c 9-`
set hash = `grep "#hash = " results.log | cut -c 9-`
set hshu = `grep "#hshu = " results.log | cut -c 9-`
set hshd = `grep "#hshd = " results.log | cut -c 9-`
set cdat = `grep "#date = " results.log | cut -c 9-`
set ctim = `grep "#time = " results.log | cut -c 9-`
set user = `grep "#user = " results.log | cut -c 9-`
set mach = `grep "#mach = " results.log | cut -c 9-`
set vers = `grep "#vers = " results.log | cut -c 9-`
set cases = `grep ${mach}_ results.log | cut -d " " -f 2 | sort -u`
set totl = `grep "#totl = " results.log | cut -c 9-`
set pass = `grep "#pass = " results.log | cut -c 9-`
set fail = `grep "#fail = " results.log | cut -c 9-`

#echo "debug ${repo}"
#echo "debug ${bran}"
#echo "debug ${hash}"
#echo "debug ${hshu}"
#echo "debug ${hshd}"
#echo "debug ${cdat}"
#echo "debug ${ctim}"
#echo "debug ${user}"
#echo "debug ${mach}"
#echo "debug ${vers}"
#echo "debug ${totl}"
#echo "debug ${pass}"
#echo "debug ${fail}"
#echo "debug ${cases}"

set xcdat = `echo $cdat | sed 's|-||g' | cut -c 3-`
set xctim = `echo $ctim | sed 's|:||g'`
set shhash = `echo $hash | cut -c 1-10`
set shrepo = `echo $repo | tr '[A-Z]' '[a-z]'`

if ("${shrepo}" !~ "*cice-consortium*") then
  set hfile = {$hfile}_forks
  set mfile = {$mfile}_forks
  set vfile = {$vfile}_forks
  set bfile = {$bfile}_forks
endif

#==============================================================
# Create results table
#==============================================================

set ofile = "${hash}.${xcdat}${xctim}"
set outfile = "${wikiname}/${tsubdir}/${ofile}.md"
mkdir -p ${wikiname}/${tsubdir}
echo "${0}: writing to ${outfile}"

if (-e ${outfile}) rm -f ${outfile}

cat >! ${outfile} << EOF

| Build | Run | Test | Regression | Compare | Timing | Case |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ |
EOF

set noglob
set green  = "\![#00C000](https://placehold.it/15/00C000/000000?text=+)"
set red    = "\![#F00000](https://placehold.it/15/F00000/000000?text=+)"
set orange = "\![#FFA500](https://placehold.it/15/FFA500/000000?text=+)"
set yellow = "\![#FFE600](https://placehold.it/15/FFE600/000000?text=+)"
set gray   = "\![#AAAAAA](https://placehold.it/15/AAAAAA/000000?text=+)"

@ ttotl = 0
@ tpass = 0
@ tfail = 0

foreach case ( ${cases} )

  @ ttotl = $ttotl + 1

  set chkpass = 1

  set fbuild = `grep " ${case} " results.log | grep " build" | cut -c 1-4`
  set fregr  = `grep " ${case} " results.log | grep " compare" | cut -c 1-4`
  set fcomp  = `grep " ${case} " results.log | grep " bfbcomp" | cut -c 1-4`
  set vregr  = `grep " ${case} " results.log | grep " compare" | cut -d " " -f 4`
  set vcomp  = `grep " ${case} " results.log | grep " bfbcomp" | cut -d " " -f 4`
  set ftime  = ""

  if (${case} =~ *_restart_*) then
    set frun   = `grep " ${case} " results.log | grep " run-initial" | cut -c 1-4`
    set frun   = `grep " ${case} " results.log | grep " run-restart" | cut -c 1-4`
    set ftest  = `grep " ${case} " results.log | grep " exact-restart" | cut -c 1-4`
  else if (${case} =~ *_smoke_*) then
    set frun   = `grep " ${case} " results.log | grep " run" | cut -c 1-4`
    set ftest  = `grep " ${case} " results.log | grep " run" | cut -c 1-4`
  endif

  set rbuild = ${yellow}
  set rrun   = ${yellow}
  set rtest  = ${yellow}
  set rregr  = ${yellow}
  set rcomp  = ${yellow}
  set rtime  = ${yellow}

  if (${fbuild} == "PASS") set rbuild = ${green}
  if (${frun}   == "PASS") set rrun   = ${green}
  if (${ftest}  == "PASS") set rtest  = ${green}
  if (${fregr}  == "PASS") set rregr  = ${green}
  if (${fcomp}  == "PASS") set rcomp  = ${green}

  if (${fbuild} == "FAIL") set rbuild = ${red}
  if (${frun}   == "FAIL") set rrun   = ${red}
  if (${ftest}  == "FAIL") set rtest  = ${red}
  if (${fregr}  == "FAIL") set rregr  = ${red}
  if (${fcomp}  == "FAIL") set rcomp  = ${red}

  if (${fbuild} == "") set rbuild = ${red}
  if (${frun}   == "") set rrun   = ${red}
  if (${ftest}  == "") set rtest  = ${red}
  if (${fregr}  == "") set rregr  = ${gray}
  if (${fcomp}  == "") set rcomp  = ${gray}
  if (${ftime}  == "") set rtime  = ${gray}

  if (${rbuild} == ${red}) set chkpass = 0
  if (${rrun}   == ${red}) set chkpass = 0
  if (${rtest}  == ${red}) set chkpass = 0

  if (${chkpass} == 1) then
     @ tpass = $tpass + 1
  else
     @ tfail = $tfail + 1
  endif

  set xcase = `echo $case | sed 's|_| |g'`
  set xvcomp = `echo $vcomp | sed 's|_| |g'`
  #echo "debug | ${rbuild} | ${rrun} | ${rtest} | ${rregr} ${vregr} | ${rcomp} ${vcomp} | ${case} |" 
  echo "| ${rbuild} | ${rrun} | ${rtest} | ${rregr} ${vregr} | ${rcomp} ${xvcomp} | ${rtime} | ${xcase} |" >> ${outfile}

end

set tcolor = ${green}
if (${tfail} > 0) set tcolor = ${yellow}
@ chk = (${ttotl} / 10)
if (${tfail} >= ${chk}) set tcolor = ${orange}
@ chk = (${ttotl} / 5)
if (${tfail} >= ${chk}) set tcolor = ${red}

unset noglob

mv ${outfile} ${outfile}.hold
cat >! ${outfile} << EOF
- repo = **${repo}** : **${bran}**
- hash = ${hash}
- hash created by ${hshu} ${hshd}
- vers = ${vers}
- tested on ${mach}, ${user}, ${cdat} ${ctim} UTC
- raw results: ${totl} total tests: ${pass} pass, ${fail} fail
- ${ttotl} total tests: ${tpass} pass, ${tfail} fail
EOF
cat ${outfile}.hold >> ${outfile}

cat >> ${outfile} << EOF

--------

EOF

#==============================================================

set hashfile = "${wikiname}/${tsubdir}/${hfile}.md"
set machfile = "${wikiname}/${tsubdir}/${mfile}.md"
set versfile = "${wikiname}/${tsubdir}/${vfile}.md"
set branfile = "${wikiname}/${tsubdir}/${bfile}.md"

foreach xfile ($hashfile $machfile $versfile $branfile)
  if (-e ${xfile}) then
    cp -f ${xfile} ${xfile}.prev
  endif
end

#=====================
# update hashfile
#=====================

set chk = 0
if (-e ${hashfile}) set chk = `grep "\*\*${hash}" ${hashfile} | wc -l`
if ($chk == 0) then
cat >! ${hashfile} << EOF
**${hash}** :

| machine | version | date/time | pass | fail | total |
| ------ | ------ | ------ | ------  | ------ | ------ |
| ${mach} | ${vers} | ${cdat} ${ctim} | ${tpass} | ${tfail} | ${tcolor} [${ttotl}](${ofile}) |

EOF
if (-e ${hashfile}.prev) cat ${hashfile}.prev >> ${hashfile}

else
  set oline = `grep -n "\*\*${hash}" ${hashfile} | head -1 | cut -d : -f 1`
  @ nline = ${oline} + 3
  sed -i "$nline a | ${mach} | ${vers} | ${cdat} ${ctim} | ${tpass} | ${tfail} | ${tcolor} [${ttotl}](${ofile}) | " ${hashfile}
endif

#=====================
# update versfile
#=====================

set chk = 0
if (-e ${versfile}) set chk = `grep "\*\*${vers}" ${versfile} | wc -l`
if ($chk == 0) then
cat >! ${versfile} << EOF
**${vers}** :

| machine | hash | date/time | pass | fail | total |
| ------ | ------ | ------ | ------  | ------ | ------ |
| ${mach} | ${shhash} | ${cdat} ${ctim} | ${tpass} | ${tfail} | ${tcolor} [${ttotl}](${ofile}) |

EOF
if (-e ${versfile}.prev) cat ${versfile}.prev >> ${versfile}

else
  set oline = `grep -n "\*\*${vers}" ${versfile} | head -1 | cut -d : -f 1`
  @ nline = ${oline} + 3
  sed -i "$nline a | ${mach} | ${shhash} | ${cdat} ${ctim} | ${tpass} | ${tfail} | ${tcolor} [${ttotl}](${ofile}) | " ${versfile}
endif

#=====================
# update machfile
#=====================

set chk = 0
if (-e ${machfile}) set chk = `grep "\*\*${mach}" ${machfile} | wc -l`
if ($chk == 0) then
cat >! ${machfile} << EOF
**${mach}** :

| version | hash | date/time | pass | fail | total |
| ------ | ------ | ------ | ------  | ------ | ------ |
| ${vers} | ${shhash} | ${cdat} ${ctim} | ${tpass} | ${tfail} | ${tcolor} [${ttotl}](${ofile}) |

EOF
if (-e ${machfile}.prev) cat ${machfile}.prev >> ${machfile}

else
  set oline = `grep -n "\*\*${mach}" ${machfile} | head -1 | cut -d : -f 1`
  @ nline = ${oline} + 3
  sed -i "$nline a | ${vers} | ${shhash} | ${cdat} ${ctim} | ${tpass} | ${tfail} | ${tcolor} [${ttotl}](${ofile}) | " ${machfile}
endif

#=====================
# update branfile
#=====================

set chk = 0
if (-e ${branfile}) set chk = `grep "\*\*${bran}" ${branfile} | wc -l`
if ($chk == 0) then
cat >! ${branfile} << EOF
**${bran}** **${repo}**:

| machine | hash | date/time | pass | fail | total |
| ------ | ------ | ------ | ------  | ------ | ------ |
| ${mach} | ${shhash} | ${cdat} ${ctim} | ${tpass} | ${tfail} | ${tcolor} [${ttotl}](${ofile}) |

EOF
if (-e ${branfile}.prev) cat ${branfile}.prev >> ${branfile}

else
  set oline = `grep -n "\*\*${bran}" ${branfile} | head -1 | cut -d : -f 1`
  @ nline = ${oline} + 3
  sed -i "$nline a | ${mach} | ${shhash} | ${cdat} ${ctim} | ${tpass} | ${tfail} | ${tcolor} [${ttotl}](${ofile}) | " ${branfile}
endif

#=====================
# update wiki
#=====================

cd ${wikiname}
git add ${tsubdir}/${ofile}.md
git add ${tsubdir}/${hfile}.md
git add ${tsubdir}/${mfile}.md
git add ${tsubdir}/${vfile}.md
git add ${tsubdir}/${bfile}.md
git commit -a -m "update $hash $mach"
git push origin master
cd ../

