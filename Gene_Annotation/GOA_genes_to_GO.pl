#!/usr/bin/perl

use strict;
use warnings;
use List::MoreUtils qw{uniq};
use Storable;
die "$0 [blast tsv file] [GOA_db]\n" unless @ARGV >= 2;
my $blast = $ARGV[0];
my $db = $ARGV[1];
my (%acc_gene, $ac_go, %gene_ac, %ac_gene);
my %go_term = (
    'MF'=>'go_function',
    'CC'=>'go_component',
    'BP'=>'go_process',
    'molecular_function'=>'go_function',
    'cellular_component'=>'go_component',
    'biological_process'=>'go_process',
    'F'=>'go_function',
    'P'=>'go_process',
    'C'=>'go_component'
);

sub get_ac_from_tsv {
	# /opt/bio/diamond blastp -d zebrafish_prot.dmnd -e 1e-5 -p 16 -q ../dj_protV2.fa --sensitive -k 1 -f 6 -o dj_to_zebrafish.tsv
	# DJ_000001-T1	tr|E7EY87|E7EY87_DANRE	100.0	30	0	0	10	39	432	461	1.9e-11	66.2
    # only the 1st and 2nd is needed.
	my $tsv = shift;
    open(my $in_fh, "$tsv") or die "$!";
    while(<$in_fh>) {
        next if /^$|^#/;
        my @p = split(/\t/);
        my $gene = (split(/\s+/, $p[0]))[0];
        my $ac   = (split(/\|/, $p[1]))[1];
        $ac_gene{$ac}   = $gene;
		$gene_ac{$gene} = $ac;
    }

}

sub mk_goa_db {
    # make goa db, will return $ac_go.
    my $db    = shift;
    my $cache = "$db.cache";
    
	if (1) {
        if (-f "$cache") {
            my $t1    = (stat($cache))[9];
            my $t2    = (stat($db))[9];
            if ( $t1 > $t2) {
                $ac_go = retrieve("$cache");
                return;
            }
        }
    }
	open(my $in_fh, "$db") or die "$!";
	while(<$in_fh>) {
		# UniProtKB	A0A023PGP0	bach1b		GO:0000976	PMID:24652768	IDA		F	BTB and CNC homology 1, basic leucine zipper transcription factor 1 b	bach1b	protein	taxon:7955	20140918	ZFIN		
		next if /^#|^!/;
		my @p  = split(/\t/);
		my $ac = $p[1];
		my $go = $p[4];
		die "No such GO = $go\n" if $go !~/GO:/;
		push @{$ac_go->{$ac}}, $go;
	}
	close $in_fh;
    store $ac_go, "$cache" if 1;
}

get_ac_from_tsv($blast);
mk_goa_db($db);

foreach my $gene (sort keys %gene_ac) {
    print $gene, "\t", join(',', uniq @{$ac_go->{$gene_ac{$gene}}}), "\n";
}
