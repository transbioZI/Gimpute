#!/bin/perl

# File   : snpsWithSamePos.pl
# Author : Dietmar Lippold
# Version: 13. September 2016

# Determines all SNPs in a first bim file which have the same position,
# i.e. the same combination of chromosome and base pair position, as a
# SNP (with the same or with another name) in a second bim file. These
# SNPs are output to STDOUT. The positions in each bim file must be unique,
# i.e. there may not be two SNPs in one of the bim files with the same
# positions.

# Usage:
# ./snpsWithSamePos.pl <bim-file-1> <bim-file-2>
#
# <bim-file-1> is the bim file from which the SNPs are determined.
# <bim-file-2> is the bim file from which the positions are compared to
# those of <bim-file-1>.


use strict;


## Config variables.


## Changeable variables

# The name of the first bim file.
my $bimFile1;

# The name of the second bim file.
my $bimFile2;

# For every position of the first bim file the name of the respective SNP
# is stored.
my %posSnps1 = ();

# For every position of the second bim file the name of the respective SNP
# is stored.
my %posSnps2 = ();


## Main    program   

if (scalar(@ARGV) != 2) {
    print STDERR "Wrong number of parameters.\n";
    print STDERR "Usage: ./snpsWithSamePos.pl <bim-file-1> <bim-file-2>\n";
    exit(-1);
}

my $bimFile1 = $ARGV[0];
my $bimFile2 = $ARGV[1];

# Read the first bim file.

open(BIMFILE, $bimFile1) or die "Can't open $bimFile1: $!";
while (<BIMFILE>) {
    # Replace DOS line endings.
    s/\r[\n]?/\n/gm;

    chomp();
    my @line = split(' ');

    my $chrNo = $line[0];
    my $bpPos = $line[3];
    my $chrPos = "${chrNo}:${bpPos}";

    if (exists($posSnps1{$chrPos})) {
        print STDERR "Warning, position is contained twice in dataset: $bimFile1, $chrPos\n";
    }

    $posSnps1{$chrPos} = $line[1];
}
close(BIMFILE);

# Read the second bim file.

open(BIMFILE, $bimFile2) or die "Can't open $bimFile2: $!";
while (<BIMFILE>) {
    # Replace DOS line endings.

    s/\r[\n]?/\n/gm;

    chomp();
    my @line = split(' ');

    my $chrNo = $line[0];
    my $bpPos = $line[3];
    my $chrPos = "${chrNo}:${bpPos}";

    if (exists($posSnps1{$chrPos})) {
        print STDERR "Warning, position is contained twice in dataset: $bimFile2, 
$chrPos\n";
    }

    $posSnps2{$chrPos} = $line[1];
}
close(BIMFILE);

# Determine and output the SNPs with the same position.

foreach my $chrPos (sort(keys(%posSnps1))) {
    if (exists($posSnps2{$chrPos})) {
        print $posSnps1{$chrPos} . "\n";
    }
}

