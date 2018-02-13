#!/bin/csh

# Check to see if test case directory was passed
if ( $1 == "" ) then
  echo "To generate timeseries plots, this script must be called with a filename."
  echo "Example: ./timeseries.csh /work/username/case01/ice_diag.itd"
  exit -1
endif

#set basename = `echo $1:t`
set basename = $1

set fieldlist=("area fraction  " \
               "avg ice thickness (m)" \
               "avg snow depth (m)" \
               "air temperature (C)" \
               "shortwave radiation sum" \
               "longwave radiation" \
               "snowfall" \
               "avg salinity (ppt)" \
               "surface temperature(C)" \
               "outward longwave flx" \
               "sensible heat flx" \
               "latent heat flx" \
#               "subl/cond (m ice)" \
               "top melt (m)" \
               "bottom melt (m)" \
               "lateral melt (m)" \
               "new ice (m)" \
               "congelation (m)" \
#               "snow-ice (m)" \
               "intnl enrgy chng(W/m^2)")

# Get the filename for the latest log
set logfile = $1

# Loop through each field and create the plot
foreach field ($fieldlist:q)
  set fieldname = `echo "$field" | sed -e 's/([^()]*)//g'`
  set search = "'$fieldname'\|istep1"
  rm -f data.txt
  # Create the new data file that houses the timeseries data
  # assumes daily output
  foreach line ("`egrep $search $logfile`")
    if ("$line" =~ *"istep1"*) then
      set argv = ( $line )
      set date = $4
      @ hour = ( $6 / 3600 )
    else
      set data1 = `echo $line | rev | cut -d ' ' -f1 | rev`
      echo "$date-$hour,$data1" >> data.txt
    endif
  end
  set format = "%Y%m%d-%H"

  # Set x-axis limits
    # User-defined x-axis limits
  # set xrange = "set xrange ['20150301':'20150901']"
    # ...Or let gnuplot determine x-axis limits
  set xrange = ""

  # Set y-axis limits
  if ("$fieldname" =~ *"area fraction"*) then
    set yrange = "set yrange [0:1]"
  else if ("$fieldname" =~ *"avg ice thickness"*) then
    set yrange = "set yrange [0:5]"  # in meters
  else if ("$fieldname" =~ *"avg snow depth"*) then
    set yrange = "set yrange [0:0.5]"  # in meters
  else
    set yrange = ""
  endif

  set output = `echo $fieldname | sed 's/ /_/g'`
  set casename = `echo $basename | rev | cut -d / -f 1-2 | rev | sed 's/\//, /'`
  set fname_base = "${casename}_${output}"
  set output_fname = "${basename}_${output}.png"

  echo "Plotting data for '$fieldname' and saving to $output_fname"

# Call the plotting routine, which uses the data in the data.txt file
gnuplot << EOF > $output_fname
# Plot style
set style data points

set datafile separator ","

# Term type and background color, canvas size
set terminal png size 1920,960

# x-axis 
set xdata time
set timefmt "$format"
set format x "%Y/%m/%d"

# Axis tick marks
set xtics rotate

set title "$fname_base"
set ylabel "$field"
set xlabel "Simulation Day"

# Set y-axis limits
$yrange

# Set x-axis limits
$xrange

# Since only 1 field is plotted, turn off legend
set key off

plot "data.txt" using (timecolumn(1)):2 with lines lw 2 lt 1 title " "

EOF

# Delete the data file
/bin/rm data.txt
end
