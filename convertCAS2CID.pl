#! /usr/bin/perl -w
###illumina ID is located at the column 14. need to unique the data 
use strict;
#my $usage = "usage: $0 <csv> \n";
#die $usage unless @ARGV == 1;
my ($filename, $col) = @ARGV;
my %data;
my %getIt;
my %result;
open (A, "<$filename")||die "could not open $filename\n";
my $tt=<A>;
for($tt){s/\r//gi;}
while (my $xx= <A>){
    chomp $xx;
    for($xx){s/\r//gi;}
    my @tmp=split /\t/, $xx;
    $data{lc($tmp[$col-1])}=$xx;
    $getIt{lc($tmp[$col-1])}=0;
    
}
close A;

print "CID\tMatchName\t$tt\n";
#open (B, "<CID-Synonym-unfiltered")||die "could not open CID-Synonym-unfiltered\n";
open (B, "<CID-CAS.1")||die "could not open CID-CAS.1\n";
while (my $x=<B>){
    chomp $x;
    for ($x){s/\r//gi;}
    my @tmp=split /\t/, $x;
    if($tmp[1] =~ /CAS-(.+)$/){
	my $id=$1;
	#print STDERR $id, "\n";
      if(defined $data{$id}){
	print $x, "\t", $tmp[1], "\t", $data{$id}, "\n";
	$getIt{$id}=1;
	$result{$id}{$tmp[0]}=1;

      }
    }
}
close B;

foreach my $i (keys %getIt){
    if($getIt{$i}==0){
	print "NA\t$i\t", $data{$i}, "\n";
    }
}

open (A, "<$filename")||die "could not open $filename\n";
open (B, ">$filename.$col.convert")||die "could not open $filename.$col.convert\n";
my $tt=<A>;
chomp $tt;
for($tt){s/\r//gi;}
print B $tt, "\tCID\n";
while (my $xx= <A>){
    chomp $xx;
    for($xx){s/\r//gi;}
    my @tmp=split /\t/, $xx;
    print B $xx, "\t";
    foreach my $i (keys %{$result{lc($tmp[$col-1])}}){
	print B   $i, ";";
    }
    print B "\n";
}
close A;
close B;
