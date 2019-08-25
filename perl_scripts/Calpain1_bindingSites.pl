#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;

#cript to find if calpain binds to Dyrk1a seq or any other sequence
################################################################################
#9c
print "####Does calpain binds to this seq?  \n";
print "###############\n\n";
#######################################	
open IN, "/Users/marmmoreno/Desktop/190218/calpain-Dyrk1a/finding_Calping_binding/Dyrk1a_mo" or die $!;
open OUT,">","/Users/marmmoreno/Desktop/190218/calpain-Dyrk1a/finding_Calping_binding/Dyrk1a_mo_calpainSites.txt" or die $!;

my $header= <IN>; #the header is saved
my $calpainCount=0;	
my $calpainCount2=0;
my $calpainCount3=0;	
my $calpainCount4=0;	
my $calpainCount5=0;	
my $calpainCount6=0;
my $resPos;
my @resSitesPos;
print OUT $header,"\n";
while(<IN>) {
my $string = $_;
while ($string=~ /([WP][LTV][KYRTFM][STRA][PS][PW]P)/g){
	$resPos= pos($string)-3; #pos give us the length until the end of the pattern found! 
	push @resSitesPos, $resPos;
	print OUT $resPos," CalpainI FULL Jin 2015 motif ", "\n";
	++ $calpainCount;
}

while ($string=~ /([WP][LTV][KYRTFM])/g){
	$resPos= pos($string);  
	push @resSitesPos, $resPos;
	print OUT $resPos," CalpainI P3-P2-P1 motif ", "\n";
	++ $calpainCount2;
}


while ($string=~ /([LTV][KYRTFM][STRA])/g){
	$resPos= pos($string);  
	push @resSitesPos, $resPos;
	print OUT  $resPos," CalpainI P2-P1-P1' motif ", "\n";
	++ $calpainCount3;
}

while ($string=~ /([WP][LTV][KYRTFM][STRA])/g){
	$resPos= pos($string);  
	push @resSitesPos, $resPos;
	print OUT $resPos," CalpainI P3-P2-P1-P1' motif ", "\n";
	++ $calpainCount4;
}
while ($string=~ /([KYRTFM][STRA])/g){
	$resPos= pos($string);  
	push @resSitesPos, $resPos;
	print OUT $resPos," CalpainI P1-P1' motif ", "\n";
	++ $calpainCount5;
}
while ($string=~ /([LTV][KYRTFM])/g){
	$resPos= pos($string);  
	push @resSitesPos, $resPos;
	print OUT  $resPos," CalpainI P2-P1 motif", "\n";
	++ $calpainCount6;
}
}
print OUT "\n";
#print "@resSitesPos\n";

close IN or die $!;
open IN, "/Users/marmmoreno/Desktop/190218/calpain-Dyrk1a/finding_Calping_binding/Dyrk1a_t588n_mo" or die $!;

$header= <IN>; 
$calpainCount=0;	
$calpainCount2=0;
$calpainCount3=0;	
$calpainCount4=0;	
$calpainCount5=0;	
$calpainCount6=0;
my $resPos2;
my @resSitesPos2;
print OUT $header,"\n";
while(<IN>) {
my $string = $_;
while ($string=~ /([WP][LTV][KYRTFM][STRA][PS][PW]P)/g){
	$resPos2= pos($string)-3; #pos give us the length until the end of the pattern found! 
	push @resSitesPos2, $resPos2;
	print OUT $resPos2," CalpainI FULL Jin 2015 motif ", "\n";
	++ $calpainCount;
}

while ($string=~ /([WP][LTV][KYRTFM])/g){
	$resPos2= pos($string); #pos give us the length until the end of the pattern found! 
	push @resSitesPos2, $resPos2;
	print OUT $resPos2," CalpainI P3-P2-P1 motif ", "\n";
	++ $calpainCount2;
}


while ($string=~ /([LTV][KYRTFM][STRA])/g){
	$resPos2= pos($string); #pos give us the length until the end of the pattern found! 
	push @resSitesPos2, $resPos2;
	print OUT  $resPos2," CalpainI P2-P1-P1' motif ", "\n";
	++ $calpainCount3;
}


while ($string=~ /([WP][LTV][KYRTFM][STRA])/g){
	$resPos2= pos($string); #pos give us the length until the end of the pattern found! 
	push @resSitesPos2, $resPos2;
	print OUT $resPos2," CalpainI P3-P2-P1-P1' motif ", "\n";
	++ $calpainCount4;
}
while ($string=~ /([KYRTFM][STRA])/g){
	$resPos2= pos($string); #pos give us the length until the end of the pattern found! 
	push @resSitesPos2, $resPos2;
	print OUT $resPos2," CalpainI P1-P1' motif ", "\n";
	++ $calpainCount5;
}
while ($string=~ /([LTV][KYRTFM])/g){
	$resPos2= pos($string); #pos give us the length until the end of the pattern found! 
	push @resSitesPos2, $resPos2;
	print OUT  $resPos2," CalpainI P2-P1 motif", "\n";
	++ $calpainCount6;
}
}

close IN or die $!;
close OUT or die $!;


#NOTE:
# Is pos  -3 instead of -4 because the system is 0 based 
# and the position that we are interested in reporting is always the -4 
# in 1 based system (thats how humans usually count! :) )


print "\n\n\n";
print "########################\n\n";
