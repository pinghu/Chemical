#! /usr/bin/perl -w
###illumina ID is located at the column 14. need to unique the data 
use strict;
#my $usage = "usage: $0 <csv> \n";
#die $usage unless @ARGV == 1;
my ($filename, $col) = @ARGV;
my %getIt;
my %data;
open (A, "<$filename")||die "could not open $filename\n";
my $tt=<A>;
chomp $tt;
for($tt){s/\r//gi;}
while (my $xx= <A>){
    chomp $xx;
    for($xx){s/\r//gi;}
    my @tmp=split /\t/, $xx;
    if(defined $tmp[$col-1]){
	my @ttt=split /\;/, $tmp[$col-1];
	for (my $i=0; $i< scalar @ttt; $i++){
	    $getIt{$ttt[$i]}{$xx}=0;
	    #print STDERR $ttt[$i], "\n";
	}
    }else{
#	print STDERR $xx, "\n";
	 $getIt{"NA"}{$xx}=0;
    }
}
close A;

print "$tt\tCASNo\tMatchCID\n";
#open (B, "<CID-CAS-07-11-2023")||die "could not open CID-CAS-07-11-2023\n";
open (B, "<CID-CAS.1")||die "could not open CID-CAS.1\n";
while (my $x=<B>){
    chomp $x;
    for ($x){s/\r//gi;}
    my @tmp=split /\t/, $x;
   
    if(defined $getIt{$tmp[0]}){
	foreach my $j (sort keys %{$getIt{$tmp[0]}}){
	    
	    print $j, "\t", $tmp[1],"\t", $tmp[0], "\n";
	    $getIt{$tmp[0]}{$j}=1;
	    $data{$j}=1;
	}
    }
    
}
close B;

foreach my $i (keys %getIt){
    foreach my $j (keys %{$getIt{$i}}){
	if($getIt{$i}{$j}==0){
	    if(! defined $data{$j}){
		print "$j\tNA\tNA\n";
		$data{$j}=1;
	    }
	}
    }
}
