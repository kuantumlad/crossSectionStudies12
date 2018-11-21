#!/bin/tcsh 

# define script parameters here
set beam_energy=$1
set num_lund_files=$2
set num_events_per_file=$3
set output_lund_dir="/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/elastic/lund"

# check if input file exists
if (-f "input$beam_energy.inp") then
    #change number of events generated per lund file
    echo " Changing number of events per file to " $num_events_per_file
    rm test_input$beam_energy.inp
    awk -v ev=$num_events_per_file 'NR==10 {$0=ev} { print }' input$beam_energy.inp > test_input$beam_energy.inp
else
    echo "File input$beam_energy.inp does not exist"
endif
    
set fNumber = 0
 while ( $fNumber <= $num_lund_files) 
 echo " creating lund file " $fNumber
 elast_gen < test_input$beam_energy.inp
 mv clas12.lund $output_lund_dir/clas12_${beam_energy}GeV_${fNumber}.lund
@ fNumber++
 end
