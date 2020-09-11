#!/usr/bin/perl
# will make gene_to_GO file.
# gene_id<tab>go1,go2...
use strict;
use warnings;
use List::MoreUtils qw{uniq};

my $gff = $ARGV[0];
open(my $in_fh, "$gff") or die "$!";
while(<$in_fh>) {
    if (/Ontology_term=/) {
        /ID=([^;]+);.+Ontology_term=([^;]+);/;
        my $id = $1;
        my $go = $2;
        my @go = split(/,/, $go);
        print join("\t", $id, join(',', @go)), "\n";
    }
}
close $in_fh;