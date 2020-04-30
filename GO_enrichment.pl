#!/usr/bin/perl
# do fisher's one-tailed exact test for overrepresent enrichment.
# from gff3 file.
use strict;
use warnings;
use Text::NSP::Measures::2D::Fisher::right;
use Statistics::Multtest qw(:all);
use Data::Dumper;
use HTTP::Tiny;
use JSON;
die "$0 [back gff] [anno gff] [GO db]\n" unless @ARGV == 3;
my $http = HTTP::Tiny->new();
my $gff = $ARGV[0];
my $ann = $ARGV[1];
my $get = $ARGV[2]; # get annotation form EMBL.
my ($full, $anno, $ref_num, $list_num, $go_db, $alt);

open(my $in_fh, "$gff") or die "$!";
while(<$in_fh>) {
    if (/Ontology_term=/) {
        /ID=([^;]+);.+Ontology_term=([^;]+);/;
        my $id = $1;
        my $go = $2;
        my @go = split(/,/, $go);
        map {$full->{$_}->{$id}++;} @go;
        $ref_num++;
    }
}
close $in_fh;

open($in_fh, "$ann") or die "$!";
while(<$in_fh>) {
    if (/Ontology_term=/) {
        /ID=([^;]+);.+Ontology_term=([^;]+);/;
        my $id = $1;
        my $go = $2;
        my @go = split(/,/, $go);
        map {$anno->{$_}->{$id}++;} @go;
        $list_num++;
    }
}
close $in_fh;
my (@out, @pval);
foreach my $go (keys %$anno) {
    #
    #           Gene_list Genome
    #in anno       a        b
    #not in anno   c        d
    #
    my @tmp_out;
    my @genes_in_list = (keys %{$anno->{$go}});
    my @genes_in_ref  = (keys %{$full->{$go}});
    my $a = scalar(@genes_in_list);
    my $b = scalar(@genes_in_ref) - $a;
    my $c = $list_num - $a;
    my $d = $ref_num  - $b - $list_num;
    my $n = $a + $b + $c + $d;
    my $p = calculateStatistic( n11=>$a,
                                n1p=>$a+$b,
                                np1=>$a+$c,
                                npp=>$n);
    push @pval, $p > 1 ? 1 : $p;
    my $go_desc  = '';
    my $go_cat   = '';
    my $go_txt   = '';
    
    if ($get eq 'embl') {
        # get from EMBL, slow!
        my $requestURL = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/'.$go;
        my $response = $http->get($requestURL, {headers => { 'Accept' => 'application/json' }});
        if ($response->{success}) {
            my $go_sca = decode_json($response->{content});
            $go_desc = $go_sca->{'results'}->[0]->{'name'};
            $go_cat  = $go_sca->{'results'}->[0]->{'aspect'};
            $go_txt  = $go_sca->{'results'}->[0]->{'definition'}->{'text'};
        }
    } elsif (-f $get) {
        if (not defined $go_db) {
            open($in_fh, "$get") or die "$!";
            while(<$in_fh>) {
                chomp;
                #go_id	go_id	Term	Ontology	Definition	Synonym	Secondary
                my @parts = split(/\t/);
                $go_db->{$parts[0]} = \@parts; # multiple relationship.
                push @{$alt->{$parts[0]}}, \@parts;
            }
        }
        warn "$go not in $go_db\n" if not defined $go_db->{$go};
        $go_desc = $go_db->{$go}->[2];
        $go_cat  = $go_db->{$go}->[3];
        $go_txt  = $go_db->{$go}->[4];
        push @out, [$go, $p, 'fdr', $go_cat, $go_desc, $a, $a+$b, join(",", @genes_in_list), $go_txt];
    }
}

# FDR correction.
my $res = BH(\@pval); my $i=0;
print join("\t", 'GO', 'Pvalue', 'FDR', 'Ontology', 'Term', 'In list', 'In background', 'genes', 'Definition'), "\n";
foreach my $line (@out) {
    $line->[2] = $res->[$i];
    print join("\t", @$line), "\n";
    $i++;
}


