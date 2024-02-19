source inputs/parameters.tcl
set output [open "carbon.dat" "w"]

set mol [molinfo top]
set nframes [molinfo $mol get numframes]

#set list_carbon {"76" "77" "79" "80" "81" "84" "85" "145" "148" "149" "152" "153" "155" "156" "320" "321" "322" "323" "326" "328" "422" "423" "425" "426" "427" "429" "476" "479" "484"}
#set list_carbon_names {}

set n_tm [llength $list_carbon]
set n_tm_ck [llength $list_carbon_names]

for {set tmi 0} {$tmi < $n_tm} {incr tmi} {

	puts "status: $tmi / $n_tm"
	set lst_tmj {}
	for {set tmj 0} {$tmj < $n_tm} {incr tmj} {

		set group1 [atomselect $mol "name CA and (resid [lindex $list_carbon $tmi])"]
		set group2 [atomselect $mol "name CA and (resid [lindex $list_carbon $tmj])"]

		set lst {}
		for {set iframe 1} {$iframe < $nframes} {incr iframe} {
			$group1 frame $iframe
			$group2 frame $iframe
			$group1 update
			$group2 update
			lappend lst "[veclength [vecsub [measure center $group1] [measure center $group2]]]"
			}

		lappend lst_tmj [expr {double([tcl::mathop::+ {*}$lst])  / [llength $lst] }]
		unset group1 
		unset group2
		unset lst
		}

	puts $output "$lst_tmj"
	unset lst_tmj
	}

close $output
exit
