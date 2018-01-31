#!/bin/perl

# File   : dataSetImpRefCompareName.pl
# Author : Dietmar Lippold
# Version: 05. August 2016

# Compares the values of SNPs in a dataset with those in the imputation
# reference files and outputs the SNPs for which there are different
# values. In this script a SNP is defined by its rs-ID name. Additionally
# values of a SNP which can be compared are the position and the alleles.

# Usage:
# ./dataSetImpRefCompareName.pl <bim-file> <imp-ref-dir> [name|pos|alleles]
#
# <bim-file> is the bim file of the dataset, <imp-ref-dir> is the directory
# with the files *.legend.gz of the imputation reference. The last parameter
# defines the test which shall be executed. For "name" the existence of the
# name of the SNPs of the dataset in the imputation reference is tested. For
# "pos" the equality of the positions of the SNPs is tested. And for
# "alleles" the existence of the alleles of the SNPs is tested.


use strict;


## Config variables.

# The pattern of the file names with the alleles of the imputation
# reference.
my $ImpFilesPattern = "*legend.gz";

# The text with which the names of SNPs begin.
my $SnpNameBegin = "rs";


## Changeable variables

# The name of the bim file.
my $bimFile;

# The name of the directory with the files *.legend.gz of the imputation
# reference.
my $impRefDir;

# The test which shall be executed.
my $test;

# The position for every SNP in the bim file.
my %bimFilePos = ();

# The first allele for every SNP in the bim file.
my %bimFileAlleles1 = ();

# The second allele for every SNP in the bim file.
my %bimFileAlleles2 = ();

# The first allele for every SNP in the imputation reference files.
my %impFilesAlleles1 = ();

# The second allele for every SNP in the imputation reference files.
my %impFilesAlleles2 = ();

# The position for every SNP in the imputation reference files.
my %impFilesPos = ();


## Main program

if (scalar(@ARGV) != 3) {
    print STDERR "Wrong number of parameters.\n";
    print STDERR "Usage: ./dataSetImpRefCompareName.pl <bim-file> <imp-ref-dir> [name|pos|alleles]\n";
    exit(-1);
}

my $bimFile = $ARGV[0];
my $impRefDir = $ARGV[1];
my $test = $ARGV[2];

if (($test ne "name") && ($test ne "pos") && ($test ne "alleles")) {
    print STDERR "Wrong value for test: $test\n";
    print STDERR "Usage: ./dataSetImpRefCompareName.pl <bim-file> <imp-ref-dir> [name|pos|alleles]\n";
    exit(-1);
}

# Read the file with the alleles of the dataset.

open(BIMFILE, $bimFile) or die "Can't open $bimFile: $!";
while (<BIMFILE>) {
    # Replace DOS line endings.
    s/\r[\n]?/\n/gm;

    chomp();
    my @line = split(' ');

    my $snpName = $line[1];
    $bimFilePos{$snpName} = $line[3];
    $bimFileAlleles1{$snpName} = $line[4];
    $bimFileAlleles2{$snpName} = $line[5];
}
close(BIMFILE);

# Read the files with the alleles of the imputation reference.

open(IMPFILES, "/bin/zcat ${impRefDir}/${ImpFilesPattern} |" ) or die "Can't open imputation reference files: $!";
while (<IMPFILES>) {
    if (/^${SnpNameBegin}/) {
        # Replace DOS line endings.
        s/\r[\n]?/\n/gm;

        chomp();
        my @line = split(' ');

        my $snpName = $line[0];
        $impFilesPos{$snpName} = $line[1];
        $impFilesAlleles1{$snpName} = $line[2];
        $impFilesAlleles2{$snpName} = $line[3];
    }
}
close(IMPFILES);

# Output the names of SNPs from the dataset which do not exist in the
# imputation reference files.

if ($test eq "name") {
    foreach my $snpName (keys(%bimFilePos)) {
        if (! exists($impFilesPos{$snpName})) {

            print "$snpName\n";
        }
    }
}

# Output the names of SNPs from the dataset which exist in the imputation
# reference files but has a different position in the imputation reference
# files.

if ($test eq "pos") {
    foreach my $snpName (keys(%bimFilePos)) {
        if (exists($impFilesAlleles1{$snpName})
            && ($bimFilePos{$snpName} != $impFilesPos{$snpName}) ) {

            print "$snpName\n";
        }
    }
}

# Output the names of SNPs from the dataset which exist in the imputation
# reference files but have different alleles there.

if ($test eq "alleles") {
    foreach my $snpName (keys(%bimFilePos)) {
        if (exists($impFilesAlleles1{$snpName})
            && exists($impFilesAlleles2{$snpName}) ) {

            if ((($bimFileAlleles1{$snpName} ne "0")
                 && ($bimFileAlleles1{$snpName} ne $impFilesAlleles1{$snpName})
                 && ($bimFileAlleles1{$snpName} ne $impFilesAlleles2{$snpName}) )
                || (($bimFileAlleles2{$snpName} ne "0")
                    && ($bimFileAlleles2{$snpName} ne $impFilesAlleles1{$snpName})
                    && ($bimFileAlleles2{$snpName} ne $impFilesAlleles2{$snpName}) )) {

                print "$snpName\n";
            }
        }
    }
}

