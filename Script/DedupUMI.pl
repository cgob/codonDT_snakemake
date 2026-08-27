#!/usr/bin/perl
###########################################################################
## Adapter trimming + UMI-based PCR deduplication.
##
## Read layout expected:
##   5'-[UMI_left][ insert ][UMI_right][3' adapter]-3'
##
## cutadapt removes the 3' adapter, then reads are deduplicated on the full
## remaining sequence (insert + both UMIs): two reads sharing the same insert
## AND the same UMI pair are PCR duplicates, so only the first is kept.
## The UMIs are stripped from the surviving reads (sequence and quality).
##
## Usage:
##   perl DedupUMI.pl <in.fastq> <adapter> <umi_left> <umi_right> <min_insert> <out.fastq>
###########################################################################

use strict;
use warnings;

my ($in, $adapter, $umi_l, $umi_r, $min_insert, $out) = @ARGV;
die "Usage: $0 <in.fastq> <adapter> <umi_left> <umi_right> <min_insert> <out.fastq>\n"
    unless defined $out;

die "umi_left must be a non-negative integer, got '$umi_l'\n"  unless $umi_l  =~ /^\d+$/;
die "umi_right must be a non-negative integer, got '$umi_r'\n" unless $umi_r  =~ /^\d+$/;
die "min_insert must be a non-negative integer, got '$min_insert'\n" unless $min_insert =~ /^\d+$/;

# A read must carry both UMIs plus a usable insert to be worth keeping.
my $min_len = $umi_l + $umi_r + $min_insert;

my $trimmed = "$out.cutadapt.tmp";

# ---- 1. adapter trimming -------------------------------------------------
my @cmd = ('cutadapt', '-a', $adapter, '-m', $min_len, '-o', $trimmed, $in);
system(@cmd) == 0 or die "cutadapt failed (exit $?): @cmd\n";

# ---- 2. dedup on (insert + UMI pair), then strip UMIs --------------------
open my $FH,  '<', $trimmed or die "Could not open $trimmed: $!\n";
open my $OUT, '>', $out     or die "Could not open $out: $!\n";

my %seen;
my ($total, $kept, $dup, $short) = (0, 0, 0, 0);

while (my $name = <$FH>) {
    my $seq  = <$FH>;
    my $plus = <$FH>;
    my $qual = <$FH>;
    die "Truncated FASTQ record at read $total in $trimmed\n" unless defined $qual;

    chomp($name, $seq, $plus, $qual);
    $total++;

    # cutadapt -m already enforces this, but a stray short read would make the
    # substr below silently return an empty or negative-length string.
    if (length($seq) < $min_len) { $short++; next; }

    # Key on the untrimmed sequence: identical insert + identical UMIs = PCR duplicate.
    if ($seen{$seq}++) { $dup++; next; }

    print $OUT $name, "\n",
               substr($seq,  $umi_l, length($seq)  - $umi_l - $umi_r), "\n",
               $plus, "\n",
               substr($qual, $umi_l, length($qual) - $umi_l - $umi_r), "\n";
    $kept++;
}

close $FH;
close $OUT;
unlink $trimmed;

printf STDERR "%s: %d reads after adapter trimming, %d kept, %d PCR duplicates, %d too short\n",
    $in, $total, $kept, $dup, $short;
