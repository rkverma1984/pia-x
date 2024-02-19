#package require math
#package require math::statistics
source vmdscripts/statistics.tcl

source inputs/parameters.tcl

set output [open "angle.dat" "w"]

set mol [molinfo top]
set nframes [molinfo $mol get numframes]

# set list_tms {"58 to 64" "65 to 74" "75 to 82" "83 to 92" "93 to 95" "96 to 101" "102 to 111" "112 to 124" "125 to 135" "136 to 151" "152 to 156" "157 to 172" "173 to 237" "238 to 246" "247 to 255" "256 to 261" "262 to 266" "267 to 274" "275 to 284" "285 to 307" "308 to 316" "317 to 328" "329 to 335" "336 to 342" "343 to 352" "353 to 359" "360 to 374" "375 to 386" "387 to 404" "405 to 417" "418 to 426" "427 to 436" "437 to 441" "445 to 454" "455 to 465" "466 to 469" "470 to 478" "479 to 483" "484 to 496" "497 to 518" "519 to 529" "530 to 540" "541 to 557" "558 to 572" "573 to 584" "585 to 600"}
# set list_tms_names {"NT" "TM1i" "TM1m" "TM1e" "TMel1" "TM2e" "TM2m" "TM2i" "TMil1" "TM3i" "TM3m" "TM3e" "TMel2" "TM4e" "TM4i" "TMil2" "TM5i" "TM5m" "TM5e" "TMel3" "TM6e" "TM6m" "TM6i" "TMil3" "TM7i" "TM7m" "TM7e" "TMel4a" "TMel4b" "TM8e" "TM8m" "TM8i" "TMil4" "TM9i" "TM9e" "TMel5" "TM10e" "TM10m" "TM10i" "TMil5" "TM11i" "TM11e" "TMel6" "TM12e" "TM12i" "CT"}


## list length
set n_tm [llength $list_tms]
set n_tm_ck [llength $list_tms_names]

proc angle { a b } {
	# Angle between two vectors
	return [expr 57.29578 * acos( [vecdot $a $b] / ([veclength $a] * [veclength $b]))]
	}


proc avg {numbers} {
   set sum 0
   foreach number $numbers {
      set sum  [expr $sum + $number]
   }
   set average [expr $sum/[llength $numbers]]
   return $average
}

# set cutoff 1.5
set cutoff 2.5

for {set tmi 0} {$tmi < $n_tm} {incr tmi} {

	puts "status: $tmi / $n_tm"

	set lst_tmj {}
	for {set tmj 0} {$tmj < $n_tm} {incr tmj} {
	#puts "TMi: $tmi; TMj: $tmj"

		if { $tmi == $tmj } {

			set FFFFF 0
			lappend lst_tmj $FFFFF

		} else {

			set groupiO [atomselect top "name O  and (resid [lindex $list_tms $tmi])"]
			set groupiC [atomselect top "name C  and (resid [lindex $list_tms $tmi])"]
			set groupjO [atomselect top "name O  and (resid [lindex $list_tms $tmj])"]
			set groupjC [atomselect top "name C  and (resid [lindex $list_tms $tmj])"]

			set lst {}
			for {set iframe 1} {$iframe < $nframes} {incr iframe} {
			
				$groupiO frame $iframe
				$groupiC frame $iframe
				$groupjO frame $iframe
				$groupjC frame $iframe

				$groupiO update
				$groupiC update
				$groupjO update
				$groupjC update

				set lst_A {}
				set lst_B {}

				for { set isel 0 } { $isel < [llength  [$groupiO get {x y z}]] } { incr isel } {
					set vevtorA [vecsub [lindex [$groupiO get {x y z}] $isel] [lindex [$groupiC get {x y z}] $isel]]
					set AA [angle $vevtorA {0 0 1}]
					lappend lst_A $AA
					}
				set muA [::math::statistics::mean $lst_A]
				set sigmaA [::math::statistics::stdev $lst_A]
				#puts "muA: $muA; sigmaA: $sigmaA"

				set final_VA {0 0 0}
				for { set isel 0 } { $isel < [llength  [$groupiO get {x y z}]] } { incr isel } {
					set AAA [angle [vecsub [lindex [$groupiO get {x y z}] $isel] [lindex [$groupiC get {x y z}] $isel]] {0 0 1}]
					if { abs([expr $AAA - $muA]) < [expr $sigmaA*$cutoff] } {
						#puts "[expr abs([expr $AAA - $muA])], [expr $sigmaA*1.5]"
						set final_VA [vecadd $final_VA [vecsub [lindex [$groupiO get {x y z}] $isel] [lindex [$groupiC get {x y z}] $isel]]]
						#puts "$final_VA"
						#lappend lst_B $AAA 
						}
					}



				for { set jsel 0 } { $jsel < [llength  [$groupjO get {x y z}]] } { incr jsel } {
					set vevtorB [vecsub [lindex [$groupjO get {x y z}] $jsel] [lindex [$groupjC get {x y z}] $jsel]]
					set BB [angle $vevtorB {0 0 1}]
					lappend lst_B $BB
					}
				set muB [::math::statistics::mean $lst_B]
				set sigmaB [::math::statistics::stdev $lst_B]
				#puts "muB: $muB; sigmaB: $sigmaB"

				set final_VB {0 0 0}
				for { set jsel 0 } { $jsel < [llength  [$groupjO get {x y z}]] } { incr jsel } {
					set BBB [angle [vecsub [lindex [$groupjO get {x y z}] $jsel] [lindex [$groupjC get {x y z}] $jsel]] {0 0 1}]
					if { abs([expr $BBB - $muB]) < [expr $sigmaB*$cutoff] } {
						#puts "[expr abs([expr $BBB - $muB])], [expr $sigmaB*1.5]"
						set final_VB [vecadd $final_VB [vecsub [lindex [$groupjO get {x y z}] $jsel] [lindex [$groupjC get {x y z}] $jsel]]]
						#puts "$final_VB"
						#lappend lst_B $BBB 
						}
					}

				set FFFFF [angle $final_VA $final_VB]
				#puts $FFFFF


				}
				## end iframe
				lappend lst_tmj $FFFFF

			}
			## end if else ...

		}
		## end tmj

	puts $output "$lst_tmj"
	unset lst_tmj
	}
	## end tmi


close $output
exit

