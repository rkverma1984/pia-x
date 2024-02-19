source vmdscripts/statistics.tcl
source inputs/parameters.tcl
set output [open "chi1.dat" "w"]

set mol [molinfo top]
set nframes [molinfo $mol get numframes]

# set list_carbon {"86" "89" "90" "96" "106" "107" "110" "111" "114" "115" "118" "165" "183" "188" "189" "192" "193" "196" "197" "338" "342" "345" "346" "349" "350" "352" "353" "362" "365" "366" "369" "370" "372" "373"}
# set list_carbon_names {"2.61" "2.64" "2.65" "EL1.50" "3.28" "3.29" "3.32" "3.33" "3.36" "3.37" "3.40" "4.57" "EL2.52" "5.38" "5.39" "5.42" "5.43" "5.46" "5.47" "6.44" "6.48" "6.51" "6.52" "6.55" "6.56" "6.58" "6.59" "7.32" "7.35" "7.36" "7.39" "7.40" "7.42" "7.43"}

set n_tm [llength $list_carbon]
set n_tm_ck [llength $list_carbon_names]

set lst_chi1 {}
for {set tmi 0} {$tmi < $n_tm} {incr tmi} {

	puts "status: $tmi / $n_tm"
	set group1 [atomselect $mol "(resid [lindex $list_carbon $tmi]) and (name N CA CB CG or name N CA CB OG or name N CA CB SG)"]
	if {[$group1 num] == 3 } {
		#unset group1
		set group1 [atomselect $mol "(resid [lindex $list_carbon $tmi]) and (name N CA CB CG1 or name N CA CB OG1)"]
		}

	set lst {}
	for {set iframe 1} {$iframe < $nframes} {incr iframe} {

		if {[$group1 num] == 4 } {
			set chi1 [measure dihed [$group1 list] frame $iframe]
			#puts "$tmi; [$group1 num]; [$group1 text]; chi: $chi1"
			} else {
			set chi1 [expr 0]			
			#puts "$tmi; [$group1 num]; [$group1 text]; chi: 0"
			}
		if {$chi1 < 0} {
			set chi1 [expr $chi1 + 360]
			}

		lappend lst "$chi1"
		}

	lappend lst_chi1 [::math::statistics::mean $lst]
	
	unset group1
	unset chi1
	}

puts "$lst_chi1"
puts $output "$lst_chi1"

close $output
exit
