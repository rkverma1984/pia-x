#!/usr/bin/perl -w

## mda needs different segid to distinguish chain id

## usage
## perl change.pl -pdb drd3_prm5.3.pdb | tee drd3_prm5.3_mod.pdb

open( PDB_IN , $ARGV[1] ) || die "Could not open file\n" ;
while( <PDB_IN> ) {

	if (( $_ =~ /^ATOM/ ) && ( substr($_, 21,1) eq "R" ) &&  ( substr($_, 72,2) eq "C2")) {
		$_ =~ s/ C2 /PROR/g;
	} elsif (( $_ =~ /^ATOM/ ) && ( substr($_, 21,1) eq "A" ) &&  ( substr($_, 72,2) eq "C2")) {
		$_ =~ s/ C2 /PROA/g;
	} elsif (( $_ =~ /^ATOM/ ) && ( substr($_, 21,1) eq "B" ) &&  ( substr($_, 72,2) eq "C2")) {
		$_ =~ s/ C2 /PROB/g;
	} elsif (( $_ =~ /^ATOM/ ) && ( substr($_, 21,1) eq "C" ) &&  ( substr($_, 72,2) eq "C2")) {
		$_ =~ s/ C2 /PROC/g;
	} elsif (( $_ =~ /^ATOM/ ) && ( substr($_, 21,1) eq "X" ) &&  ( substr($_, 72,2) eq "C2")) {
		$_ =~ s/ C2 /PROX/g;
	}
	print $_ ;
	

}

close( PDB_IN ) ;
