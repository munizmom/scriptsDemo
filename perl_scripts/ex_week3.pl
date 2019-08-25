#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;




#######################################	
#7a
print "####Ex 7.a: \n";
print "###############\n\n";
#######################################	
#sub log2 {
#        my $n = shift;
#        return log($n)/log(2);
#    }


open IN, "/Users/marmmoreno/Documents/phdCourses/perl/week3/raw_data.txt" or die $!;
open OUT,">","/Users/marmmoreno/Documents/phdCourses/perl/week3/filtered_dataTrial.txt" or die $!;

my $valA;
my $valB;
my $valD;
my @column;
my @myImpLines;
my @myResults;
my $number = 0; 
my $header= <IN>; #the header is saved
my $ID;
my $counts=0;
while(<IN>) {
	chomp;
  	my $line = $_;
  	@column = split( /\t/, $line );
    $ID= $column[0]; 
	$valA =(log($column[1])/log(2));
  	$valB =(log($column[2])/log(2));
  	$valD =abs(log($column[1]/$column[2]))/(log(2));
# 	if ($valA > 2 or $valB > 2 or $valD > 3){
 	if (($valA > 2 or $valB > 2) and ($valD > 3)){
		push @myImpLines,$_;
		print OUT $number, "\t",$ID,"\t",$valA,"\t",$valB, "\n";
		++$number;
		#print OUT @myImpLines, "\n";
	}

}
print OUT @myResults;
print "$number";
#print "@column\n"; #ook
#print $ID;  #ook
close IN or die $!;
close OUT or die $!;

print "\n\n";
print "###############\n\n";
#######################################	

#######################################	
#7b
print "####Ex 7.b: \n";
print "###############\n\n";
#######################################	
#defining user passwd rlationships:
my $usrName;
my $Passwd;
my $tries=1;
my %usrPsswd = (
    'laura' => 'Myhome2@',
    'lith' => 'Neuro22M',
    'anne' => 'Rudolf_christmas',
    'robert' => 'Obiwan1976',
    'david' => 'DarkSide85',
    'basil' => 'mouseDet_85#',
    'mar' => 'pinky22');

print "Introduce your username: \n";
$usrName = <STDIN>;
chomp $usrName;
$usrName= lc $usrName;
	while ($tries < 6) {
		if (exists $usrPsswd{$usrName}){
		print "\nIntroduce your pasword: \n";
		$Passwd = <STDIN>;
		chomp $Passwd;
    	if ($Passwd eq $usrPsswd{$usrName}){
    		print "You are allow to know the secret: You are ... is a secret! ;)";
    	last;
    	} else{
    		print "Introduce your CORRECT pasword: \n you can try again\n\n ";
    		++$tries;
    	}
    } else {
		print "Your user account does not exist. \n Talk with the admin to request an account.\n";
	} 
}
if($tries = 5){
	print "You are not allowed any more tries. You tried 5 times! \n"
}


print "\n\n";
#######################################	
#7c
print "####Ex 7.c: \n";
print "###############\n\n";
#######################################	

open IN, "/Users/marmmoreno/Documents/phdCourses/perl/week3/active_genes.txt" or die $!;












open OUT,">","/Users/marmmoreno/Documents/phdCourses/perl/week3/filtered_data.txt" or die $!;

open IN, "/Users/marmmoreno/Documents/phdCourses/perl/week3/enhancers.txt" or die $!;
open OUT,">","/Users/marmmoreno/Documents/phdCourses/perl/week3/filtered_data.txt" or die $!;


#######################################	
#8
print "####Ex 8: \n";
print "###############\n\n";
#run it using argv: perl week3.pl /Users/marmmoreno/Documents/phdCourses/perl/week3/
#######################################	
die "No directory was specified!\n" unless (@ARGV);
my $dirpath=$ARGV[0];
chdir $dirpath or die "Cannot open file, or it does not exits or is nto a directory $!";

my $filesNumber;
my $filesNrw;
my $subdirecs;
my @files = <*>;
my $file;

while ($file = <*>){
	if (-d $file){
		++$subdirecs;
	}
	elsif (-w $file and -r $file) {
		++$filesNrw;
	}

	}

print "The number of total files is: ", scalar @files, "\n";
print "The number of rw files is: $filesNrw\n";
print "The number of subdirectories is: $subdirecs\n";

#c
print "\n########################\n\n";
use Data::Dumper;


my @ext;
my @namefile;
my %filesTypes;
my @infoFiles;
my $file;
my $line;


while ($file= <*>){
chomp $file;  #file is already the line
	@infoFiles = split( /\./, $file ); #the dot has to be scaped!
    push @ext, @infoFiles[1];
    push @namefile, $file;    	
    }

#print "extensions\n";
#print "@ext\n";
#print "fileNames\n";
#print "@namefile\n";
for my $idx (0 .. $#ext) {
	push @{$filesTypes{$ext[$idx]}},
	$namefile[$idx];
}

print Dumper \%filesTypes;


print "\n\n\n";
print "########################\n\n";



print "Exercises third week finished! :)";

