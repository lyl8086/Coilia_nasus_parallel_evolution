#!/usr/bin/env perl
# will calculate genome coverage based the given LD extent and SNP positions.
# need bedtools in path.
use strict;
use warnings;
die "$0 [snp_index] [genome_sizes] [ld_in_bp]\n" unless @ARGV == 3;
my $snp_index = $ARGV[0];
my $size      = $ARGV[1];
my $ld        = $ARGV[2];
my (%size, $sum);

# chromsome sizes.
open(my $in_fh, "$size") or die "$!";
while(<$in_fh>) {
	chomp;
	next if /^#|^$/;
	my @p = split;
	$size{$p[0]} = $p[1];
	$sum += $p[1]; # total genome size.
}
close $in_fh;

# snp index.
open($in_fh, "$snp_index") or die "$!";
open(my $out_fh, ">$snp_index.plus_ld") or die "$!";
while(<$in_fh>) {
	chomp;
	next if /^#|^$/;
	my @p = split;
	my $min = ($p[1] - $ld) <= 0 ? 1 : ($p[1] - $ld);
	my $max = ($p[1] + $ld) >= $size{$p[0]} ? $size{$p[0]} : ($p[1] + $ld);
	print $out_fh join("\t", $p[0], $min, $max), "\n";
}
close $in_fh;

# merge overlap positions.
`bedtools merge -i $snp_index.plus_ld >$snp_index.$ld.extend`;

# calculate coverage.
my $res = `awk -v s=$sum '{c+=(\$3-\$2+1)} END {print "cover bases: "c\", cover ratio=\"c/s}' $snp_index.$ld.extend`;
print $res;
