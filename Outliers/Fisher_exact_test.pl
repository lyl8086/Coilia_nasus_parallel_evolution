#!/bin/perl
# Fisher Exact Test

#      pop1	pop2
# ref   a    b
# alt   c    d
#

#
#       word2   ~word2
# word1    n11      n12 | n1p
#~word1    n21      n22 | n2p
#          --------------
#          np1      np2   npp
#


use strict;
use warnings;
use List::MoreUtils qw{uniq};
use Text::NSP::Measures::2D::Fisher::twotailed;
die "Usage:$0 [vcf] [popmap] [reference_pop_name]\n" unless @ARGV >= 3;
my $in  = $ARGV[0];
my $pop = $ARGV[1];
my $reference = $ARGV[2];
eval { &main() }; # try catch.
if ($@) {
    
    print STDERR "Error: $@\n";
    
}

sub calc_frq {
    my ($gts, $ref, $alt) = @_;
    my $allele->{$ref} = 0;
       $allele->{$alt} = 0;
       $allele->{'het'}= 0;
    my $obs_hets = 0;
    my $obs_hom1 = 0;
    my $obs_hom2 = 0;
    
    foreach my $gt (@$gts) {
    
        # unphased or phased.
        if ($gt eq '0/0' || $gt eq '0|0') {$allele->{$ref} += 2; $obs_hom1++;}
        if ($gt eq '1/1' || $gt eq '1|1') {$allele->{$alt} += 2; $obs_hom2++;}
        if ($gt eq '0/1' || $gt eq '0|1' || $gt eq '1|0' || $gt eq '1/0') {
            $allele->{$ref}++; 
            $allele->{$alt}++;
            $allele->{'het'}++;
            $obs_hets++;
        }
    }
    my $a = $allele->{$ref};
    my $b = $allele->{$alt};
    return($a, $b);
}

sub calc_factorial {
    my $k = shift;
    my $r = 1;
    for my $i(2..$k) {$r = $r*$i;}
    return $r;
}

sub calc_p {
    my ($a,$b,$c,$d,$n) = @_;
    my $x = calc_factorial($a+$b) * calc_factorial($c+$d) * 
    calc_factorial($a+$c) * calc_factorial($b+$d);
    my $y = calc_factorial($a) * calc_factorial($b) * 
    calc_factorial($c) * calc_factorial($d) * calc_factorial($n);
    return ('NA') if $y == 0;
    return ($x/$y);
}

sub calc_maf {

    my ($gts, $ref, $alt) = @_;
    my $allele->{$ref} = 0;
       $allele->{$alt} = 0;
       $allele->{'het'}= 0;
    my $obs_hets = 0;
    my $obs_hom1 = 0;
    my $obs_hom2 = 0;
    
    foreach my $gt (@$gts) {
    
        # unphased or phased.
        if ($gt eq '0/0' || $gt eq '0|0') {$allele->{$ref} += 2; $obs_hom1++;}
        if ($gt eq '1/1' || $gt eq '1|1') {$allele->{$alt} += 2; $obs_hom2++;}
        if ($gt eq '0/1' || $gt eq '0|1' || $gt eq '1|0' || $gt eq '1/0') {
            $allele->{$ref}++; 
            $allele->{$alt}++;
            $allele->{'het'}++;
            $obs_hets++;
        }
    }
    my $n    = ($allele->{$ref} + $allele->{$alt}) / 2;
    my $p    = $allele->{$ref} / (2 * $n);
    my $flag = $p >= 0.5 ? $ref : $alt;
    my $maf  = $p >= 0.5 ? $p : (1 - $p);
    my $mac  = $p >= 0.5 ? $allele->{$ref} : $allele->{$alt};
    return(sprintf("%.5f", $maf), $mac, $n, $flag);
    
}

sub calc_maf_ref {

    my ($gts, $ref, $ori) = @_;
    my $alt = 'alt';
    my $allele->{$ref} = 0;
       $allele->{$alt} = 0;
       $allele->{'het'}= 0;
       
    foreach my $gt (@$gts) {
    
        # unphased or phased.
        
        if ($gt eq '0/0' || $gt eq '0|0') {$allele->{$ref} += 2;}
        if ($gt eq '1/1' || $gt eq '1|1') {$allele->{$alt} += 2;}
        if ($gt eq '0/1' || $gt eq '0|1' || $gt eq '1|0' || $gt eq '1/0') {
            $allele->{$ref}++; 
            $allele->{$alt}++;
            $allele->{'het'}++;
        }
    }
    my $n    = ($allele->{$ref} + $allele->{$alt}) / 2;
    my $p    = $allele->{$ref} / (2 * $n);
    my $q    = $allele->{$alt} / (2 * $n);
    
    if ($ref ne $ori) {
    # if the major allele is reverse.
        $p = $q;
        $allele->{$ref} = $allele->{$alt};
    }
    return(sprintf("%.5f", $p), $allele->{$ref}, $n, $ref);    
}

sub main {
    my (@header, @od, @order, %samples, $pops);
    my $head = 0;
    open (my $in_fh, "$in") or die;
    while(<$in_fh>) {
    #################################################################
    
    #VCF format:
    #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Indiv_genotypes ...
    #Genotypes:GT:PL:DP:SP:GQ ...
    
    #################################################################
    next if /^##|^$/;
    next if $head > 0 && /^#/; 
    
    if ($head ==0 && /^#/) {
        chomp;$head++;
        push @header, split(/\t/);
        for (my $i=9; $i<=$#header; $i++) {
            $samples{$header[$i]} = $i;
        }
        open(my $pop_fh, "$pop") or die "No PopMap file!";
        while (<$pop_fh>) {
            next if /^#|^$/;
            $_ =~ s/[\r\n]|^\s+|\s+$//g;
            my @part = split;
            my $cpop = $part[1];
            my $cind = $part[0];
            die "No such individual $cind in VCF\n" if not defined $samples{$cind};
            push @od, $cpop;
            push @{$pops->{$cpop}}, $samples{$cind}; #Pop name => @indiv rank. 
            
        }
        close $pop_fh;
        
        # reorder populations, reference first.
        map {push (@order, $_) if $_ ne $reference;} uniq(@od);
        print join("\t", '#CHROM','POS', 'REF', 'ALT', @order), "\n";
        unshift(@order, $reference);
        next;
        
    }
    
    chomp;   
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genos) = split(/\t/);
    my ($fmt, $major, @out);
    my @formats = split(/:/, $format);              # Order of format:
    map { $fmt->{$formats[$_]} = $_;} 0..$#formats; # Geno => Order.
    #Iteration each population.
    my ($a, $b, $c, $d, $n, $p);
    
    if (1) {
        my @gt = ();
        foreach my $rank(@{$pops->{$reference}}) {
            #each individual of every population.
            my $geno0= $genos[$rank-9]; # Genotype of original.
            my @geno = split(/:/, $geno0);
            if (not defined $geno[$fmt->{'GT'}]) { die "GT field is requred!\n";}
            if (not defined $geno[$fmt->{'DP'}]) { die "DP field is requred!\n";}
            my $GT   = $geno[$fmt->{'GT'}];
            # Missing is skipped.
            if ($GT eq './.' || $GT eq '.' || $GT eq '.|.') {next;}
            push (@gt, $GT);
        }
        ($a,$c) = calc_frq(\@gt, $ref, $alt);
    }
    
    foreach my $name (@order) {

        next if $name eq $reference;        
        my @gt    = ();
        
        foreach my $rank(@{$pops->{$name}}) {
            #each individual of every population.
            my $geno0= $genos[$rank-9]; # Genotype of original.
            my @geno = split(/:/, $geno0);
            if (not defined $geno[$fmt->{'GT'}]) { die "GT field is requred!\n";}
            if (not defined $geno[$fmt->{'DP'}]) { die "DP field is requred!\n";}
            my $GT   = $geno[$fmt->{'GT'}];
            # Missing is skipped.
            if ($GT eq './.' || $GT eq '.' || $GT eq '.|.') {next;}
            push (@gt, $GT);
        }
        ($b,$d) = calc_frq(\@gt, $ref, $alt);
        $n = $a + $b + $c + $d;
        $p = calculateStatistic( n11=>$a,
                                 n1p=>$a+$b,
                                 np1=>$a+$c,
                                 npp=>$n);
        push @out, $p;
    }

    print join("\t", $chrom, $pos, $ref, $alt, @out), "\n";
    
    }
}
