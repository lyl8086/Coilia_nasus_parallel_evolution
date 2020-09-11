use strict;
use warnings;
my $gff = $ARGV[0];
my %GO;
open(my $in_fh, "$gff") or die "$!";
while(<$in_fh>) {
    if (/Ontology_term=/) {
        /ID=([^;]+);.+Ontology_term=([^;]+);/;
        my $id = $1;
        my $go = $2;
        my @go = split(/,/, $go);
        map {$GO{$_}++;} @go;
    }
}
close $in_fh;

my @out = (keys %GO);

print join("\n", @out);

