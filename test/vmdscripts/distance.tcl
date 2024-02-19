source inputs/parameters.tcl
set output [open "distance.dat" "w"]

set mol [molinfo top]
set nframes [molinfo $mol get numframes]

# set list_tms {"58 to 64" "65 to 74" "75 to 82" "83 to 92" "93 to 95" "96 to 101" "102 to 111" "112 to 124" "125 to 135" "136 to 151" "152 to 156" "157 to 172" "173 to 237" "238 to 246" "247 to 255" "256 to 261" "262 to 266" "267 to 274" "275 to 284" "285 to 307" "308 to 316" "317 to 328" "329 to 335" "336 to 342" "343 to 352" "353 to 359" "360 to 374" "375 to 386" "387 to 404" "405 to 417" "418 to 426" "427 to 436" "437 to 441" "445 to 454" "455 to 465" "466 to 469" "470 to 478" "479 to 483" "484 to 496" "497 to 518" "519 to 529" "530 to 540" "541 to 557" "558 to 572" "573 to 584" "585 to 600"}
# set list_tms_names {"NT" "TM1i" "TM1m" "TM1e" "TMel1" "TM2e" "TM2m" "TM2i" "TMil1" "TM3i" "TM3m" "TM3e" "TMel2" "TM4e" "TM4i" "TMil2" "TM5i" "TM5m" "TM5e" "TMel3" "TM6e" "TM6m" "TM6i" "TMil3" "TM7i" "TM7m" "TM7e" "TMel4a" "TMel4b" "TM8e" "TM8m" "TM8i" "TMil4" "TM9i" "TM9e" "TMel5" "TM10e" "TM10m" "TM10i" "TMil5" "TM11i" "TM11e" "TMel6" "TM12e" "TM12i" "CT"}

set n_tm [llength $list_tms]
set n_tm_ck [llength $list_tms_names]

for {set tmi 0} {$tmi < $n_tm} {incr tmi} {

	puts "status: $tmi / $n_tm"
	set lst_tmj {}
	for {set tmj 0} {$tmj < $n_tm} {incr tmj} {

		set group1 [atomselect $mol "name CA and (resid [lindex $list_tms $tmi])"]
		set group2 [atomselect $mol "name CA and (resid [lindex $list_tms $tmj])"]

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
