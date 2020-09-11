#!/usr/bin/perl
# extract go from uniprot accession
#
# tsv file generated using:
# blastxmlparser \
# -n 'hit.parent.query_def,hsp.evalue,hsp.bit_score,hit.accession,hit.hit_def'
#
# 1. mapping gene_id with uniprot ac from the blast tsv file
# 2. retrive go from mysql database using uniprot ac
# 3. reformat output using GO.db.tsv

use strict;
use warnings;
use DBI;
use List::MoreUtils qw{uniq};
use Storable;
die "$0 [blast tsv file]\n" unless @ARGV >= 2;
my $blast     = $ARGV[0];
my $use_cache = $ARGV[1];
my (%acc_gene, $gene_go, $go_anno);
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
my $driver="DBI:mysql";
my $database="idmapping";
my $user="liyulong";
my $host="localhost";
my $passwd="123";
my $rules="go";
my $dbh = DBI->connect("$driver:database=$database;host=$host;user=$user;password=$passwd")
or die "Can't connect: " . DBI->errstr;

sub get_go {
    # get go from ac list.
    my $ac    = shift;
    my $cache = "$blast.cache";
    if ($use_cache) {
        if (-f "$cache") {
            my $t1    = (stat($cache))[9];
            my $t2    = (stat($blast))[9];
            if ( $t1 > $t2) {
                $gene_go = retrieve("$cache");
                return;
            }
        }
    }
    foreach my $acc (@$ac) {
        my $sql = "select $rules from uniprot_select where uniac='$acc'";
        #print STDERR $sql;
        my $sth=$dbh->prepare($sql);
        $sth->execute() or die "Can't prepare sql statement". $sth->errstr;
        my @recs   = $sth->fetchrow_array;
        #print STDERR "  $recs[0]\n";
        my $gos    = $recs[0];
        if ($gos && $gos =~ /GO:/) {
            my @p  = split(/;/, $gos);
            my @GO = ();
            map {$_=~s/^\s+|\s+$//g; push @GO, $_;} @p;
            push @{$gene_go->{$acc_gene{$acc}}}, @GO;
        } 
        #$sth->finish();
    }
    store $gene_go, "$cache" if $use_cache;
}

sub get_ac_from_tsv {
    #1       1       1       Hualu_000001-T1 Hualu_000001    0.0     1603.6  A0A385FP91      Sodium/hydrogen exchanger OS=Lateolabrax maculatus OX=315492 GN=NHE5 PE=2 SV=1
    #4       2       1       Hualu_000004-T1 Hualu_000004    1.3e-52 215.7   A0A3B4TVC7      Formin_GBD_N domain-containing protein OS=Seriola dumerili OX=41447 PE=4 SV=1
    #
    my $tsv = shift;
    open(my $in_fh, "$tsv") or die "$!";
    while(<$in_fh>) {
        next if /^$|^#/;
        my @p = split(/\t/);
        my $gene = (split(/\s+/,$p[3]))[0];
        my $ac   = $p[6];
        $acc_gene{$ac} = $gene;
    }

}

get_ac_from_tsv($blast);
my @ac = (keys %acc_gene);
get_go(\@ac);

foreach my $gene (sort keys %$gene_go) {
    print $gene, "\t", join(',', uniq @{$gene_go->{$gene}}), "\n";
}