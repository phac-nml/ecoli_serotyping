#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

my $inFile = $ARGV[0];
open(my $inFH, '<', $inFile) or die "$!";

while(my $line = $inFH->getline){
    if($line =~ m/^>/){
        $line =~ s/\R//g;
        $line =~ s/gi:/gi_/g;
        $line =~ s/gi\|/gi_/g;
        $line =~ s/\s/_/g;

        $line .= "\n";
    }
    print $line;
}