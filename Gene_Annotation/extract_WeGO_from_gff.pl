use strict;
use warnings;
my $gff = $ARGV[0];
open(my $in_fh, "$gff") or die "$!";
while(<$in_fh>) {
    if (/Ontology_term=/) {
        /ID=([^;]+);.+Ontology_term=([^;]+);/;
        my $id = $1;
        my $go = $2;
        my @go = split(/,/, $go);
        print join("\t", $id, @go), "\n";
    }
}
close $in_fh;