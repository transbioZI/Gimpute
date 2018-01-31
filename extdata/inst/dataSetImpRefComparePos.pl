#!/bin/perl

# File   : dataSetImpRefComparePos.pl
# Author : Dietmar Lippold
# Version: 05. August 2016

# Compares the values of SNPs in a dataset with those in the imputation
# reference files and outputs the SNPs for which there are different
# values. In this script a SNP is defined by its chromosome and its base
# pair position. Additionally values of a SNP which can be compared are
# the alleles and the name.

# Usage:
# ./dataSetImpRefComparePos.pl <bim-file> <imp-ref-dir> [pos|alleles|name]
#
# <bim-file> is the bim file of the dataset, <imp-ref-dir> is the directory
# with the files *.legend.gz of the imputation reference. The last parameter
# defines the test which shall be executed. For "pos" the existence of
# the a SNP (i.e. combination of chromosome and base pair position) of the
# dataset in the imputation reference is tested. For "alleles" the existence
# of the alleles the SNPs is tested. And for "name" the existence of the
# rs-name of a SNP is tested.


use strict;


## Config variables.

# The pattern of the file names with the alleles of the imputation
# reference.
my $ImpFilesPatternEnd = "legend.gz";

# For the numbers of the allosomes stores its name in the file names
# with the alleles of the imputation reference. An empty string means
# that there are no files for that number.
my %alloChrName = (23 => "_chrX_nonPAR_", 24 => "", 25 => "_chrX_PAR");

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

# The name of every SNP in the bim file.
my %bimFileSnpName = ();

# The first allele for every SNP in the bim file.
my %bimFileAlleles1 = ();

# The second allele for every SNP in the bim file.
my %bimFileAlleles2 = ();

# The name of every SNP in the imputation reference files.
my %impFilesSnpName = ();

# The first allele for every SNP in the imputation reference files.
my %impFilesAlleles1 = ();

# The second allele for every SNP in the imputation reference files.
my %impFilesAlleles2 = ();


## Main program

if (scalar(@ARGV) != 3) {
    print STDERR "Wrong number of parameters.\n";
    print STDERR "Usage: ./dataSetImpRefComparePos.pl <bim-file> <imp-ref-dir> [pos|alleles|name]\n";
    exit(-1);
}

my $bimFile = $ARGV[0];
my $impRefDir = $ARGV[1];
my $test = $ARGV[2];

if (($test ne "pos") && ($test ne "alleles") && ($test ne "name")) {
    print STDERR "Wrong value for test: $test\n";
    print STDERR "Usage: ./dataSetImpRefComparePos.pl <bim-file> <imp-ref-dir> [pos|alleles|name]\n";
    exit(-1);
}

# Read the file with the alleles of the dataset.

open(BIMFILE, $bimFile) or die "Can't open $bimFile: $!";
while (<BIMFILE>) {
    # Replace DOS line endings.
    s/\r[\n]?/\n/gm;

    chomp();
    my @line = split(' ');

    my $chrNo = $line[0];
    my $bpPos = $line[3];
    my $chrPos = "${chrNo}:${bpPos}";

    $bimFileSnpName{$chrPos} = $line[1];
    $bimFileAlleles1{$chrPos} = $line[4];
    $bimFileAlleles2{$chrPos} = $line[5];
}
close(BIMFILE);

# Read the files with the alleles of the imputation reference.

for (my $chrNo = 1; $chrNo <= 25; $chrNo++) {
    my $chrName = "";
    if ($chrNo <= 22) {
        $chrName = "_chr${chrNo}_";
    } else {
        $chrName = $alloChrName{$chrNo};
    }

    if ($chrName ne "") {
        my $impFilesPattern = "*${chrName}*${ImpFilesPatternEnd}";

        open(IMPFILES, "/bin/zcat ${impRefDir}/${impFilesPattern} |" ) or die "Can't open imputation reference files: $!";
        while (<IMPFILES>) {
            if (/^${SnpNameBegin}/) {
                # Replace DOS line endings.
                s/\r[\n]?/\n/gm;

                chomp();
                my @line = split(' ');

                my $bpPos = $line[1];
                my $chrPos = "${chrNo}:${bpPos}";

                $impFilesSnpName{$chrPos} = $line[0];
                $impFilesAlleles1{$chrPos} = $line[2];
                $impFilesAlleles2{$chrPos} = $line[3];
            }
        }
        close(IMPFILES);
    }
}

# Output the SNP names from the dataset for which the combination of
# chromosome and base pair position is not stored in the imputation
# reference files.

if ($test eq "pos") {
    foreach my $chrPos (keys(%bimFileSnpName)) {
        if (! exists($impFilesSnpName{$chrPos})) {

            my $snp = $bimFileSnpName{$chrPos};
            print "$snp\n";
        }
    }
}

# Output the SNPs from the dataset for which the combination of chromosome
# and base pair position exist in the imputation reference files but has
# different alleles there.

if ($test eq "alleles") {
    foreach my $chrPos (keys(%bimFileSnpName)) {
        if (exists($impFilesAlleles1{$chrPos})
            && exists($impFilesAlleles2{$chrPos}) ) {

            if ((($bimFileAlleles1{$chrPos} ne "0")
                 && ($bimFileAlleles1{$chrPos} ne $impFilesAlleles1{$chrPos})
                 && ($bimFileAlleles1{$chrPos} ne $impFilesAlleles2{$chrPos}) )
                || (($bimFileAlleles2{$chrPos} ne "0")
                    && ($bimFileAlleles2{$chrPos} ne $impFilesAlleles1{$chrPos})
                    && ($bimFileAlleles2{$chrPos} ne $impFilesAlleles2{$chrPos}) )) {

                my $snpName = $bimFileSnpName{$chrPos};
                print "$snpName\n";
            }
        }
    }
}

# Output the SNPs from the dataset for which the combination of chromosome
# and base pair position exist in the imputation reference files but has
# an different name there.

if ($test eq "name") {
    foreach my $chrPos (keys(%bimFileSnpName)) {
        if (exists($impFilesSnpName{$chrPos})) {
            my $bimSnpName = $bimFileSnpName{$chrPos};
            my $impSnpName = $impFilesSnpName{$chrPos};

            if ($impSnpName ne $bimSnpName) {
                print "$bimSnpName\n";
            }
        }
    }
}

