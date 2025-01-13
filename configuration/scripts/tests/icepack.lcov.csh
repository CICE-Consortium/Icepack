
#--- icepack.lcov.csh ---

echo ${lcovalist}
lcov ${lcovalist} -o total.info

set lcovrepo = apcraig.github.io
set lcovhtmldir = lcov_icepack_${report_name}
if ("${useparser}" =~ "true") then
  #local parsing script
  mkdir ./${lcovhtmldir}
  ./parse_lcov.sh >! ./${lcovhtmldir}/${lcovhtmldir}.txt
  cat >! ./${lcovhtmldir}/index.html <<EOF

<embed type="text/plain" src="${lcovhtmldir}.txt"  width="800" height="20000">

EOF
else
  #genhtml
  genhtml -o ./${lcovhtmldir} --precision 2 -t "${report_name}" total.info
endif

rm -r -f ${lcovrepo}
git clone https://github.com/apcraig/${lcovrepo}
cp -p -r ${lcovhtmldir} ${lcovrepo}/

cd ${lcovrepo}
set covp0 = `grep message coverage_icepack.json | cut -d : -f 2 | cut -d \" -f 2 | cut -d % -f 1`
if ("${useparser}" =~ "true") then
  set covp  = `grep "**TOTAL**" ${lcovhtmldir}/${lcovhtmldir}.txt | cut -d \% -f 1`
else
  set covp  = `grep -i headerCovTableEntry ${lcovhtmldir}/index.html | grep % | head -1 | cut -d \> -f 2 | cut -d % -f 1`
endif
set covpi = `echo $covp | cut -d . -f 1`

set lcovhtmlname = "${covpi}%:${report_name}"
set oline = `grep -n "add_icepack_entry_here" index.html | head -1 | cut -d : -f 1`
@ nline = ${oline}
sed -i "$nline a    <li><a href="${lcovhtmldir}/index.html">${lcovhtmlname}</a></li> " index.html

set covpcolor = red
if (${covpi} > 50) set covpcolor = orange
if (${covpi} > 60) set covpcolor = yellow
if (${covpi} > 70) set covpcolor = yellowgreen
if (${covpi} > 80) set covpcolor = green
if (${covpi} > 90) set covpcolor = brightgreen

cp coverage_icepack.json coverage_icepack.json.old
cat >! coverage_icepack.json <<EOF

  {
  "schemaVersion": 1,
  "label": "lcov",
  "message": "${covp}%",
  "color": "${covpcolor}"
  }

EOF

git add .
git commit -m "add ${lcovhtmldir}"
git push origin master

