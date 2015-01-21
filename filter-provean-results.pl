#!/usr/bin/env perl
# Mike Covington
# created: 2015-01-20
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use List::Util 'sum';

my ( $gene, $list, $range, $help, $verbose );
my ( $gene, $list, $range, $gff_file, $help, $verbose );
my $out_dir = '.';
my $options = GetOptions(
    "gene=s"     => \$gene,
    "list=s"     => \$list,
    "range=s"    => \$range,
    "gff_file=s" => \$gff_file,
    "out_dir=s"  => \$out_dir,
    "help"       => \$help,
    "verbose"    => \$verbose,
);

validate_options( $gene, $list, $range, $gff_file, $out_dir, $help );

my $gene_list;
if ( defined $gene ) {
    push @$gene_list, $gene;
}
elsif ( defined $list ) {
    $gene_list = get_genes_in_list($list);
}
elsif ( defined $range ) {
    $gene_list = get_genes_in_range( $range, $gff_file );
}

exit;

sub do_or_die {
    my ( $help, $errors ) = @_;

    if ($help) {
        die usage();
    }
    elsif (@$errors) {
        my $error_string = join "\n", map {"ERROR: $_"} @$errors;
        die usage(), $error_string, "\n\n";
    }
}

sub get_genes_in_list {
    my $gene_list_file = shift;

    open my $gene_list_fh, "<", $gene_list_file;
    my @gene_list = <$gene_list_fh>;
    close $gene_list_fh;
    chomp @gene_list;

    return \@gene_list;
}

sub get_genes_in_range {
    my ( $range, $gff_file ) = @_;
    my ( $chr, $start, $end ) = split /[:-]/, $range;

    my $genes = get_gene_models($gff_file);

    my @gene_list;
    for my $mrna ( sort keys %{ $$genes{$chr} } ) {
        my $mrna_start = $$genes{$chr}{$mrna}{cds}->[0]->{start};
        my $mrna_end   = $$genes{$chr}{$mrna}{cds}->[-1]->{end};
        push @gene_list, $mrna if $mrna_end >= $start && $mrna_start <= $end;
    }

    return \@gene_list;
}

sub usage {
    return <<EOF;

Usage: $0 [options] <--filter_method> <filter_argument>

Filter Methods (Specify only one):
  --gene           Specify a single gene (--gene Solyc08g065870.2.1)
  -l, --list       Specify a file with a list of genes (one gene per line)
  -r, --range      Specify a genomic range (--range SL2.40ch01:619126-1047293)

Options:
  --gff_file       GFF3 annotation file (required when filtering by '--range')
  -o, --out_dir    Output directory [.]
  -v, --verbose    Print filter summary to STDOUT
  -h, --help       Display this usage information

EOF
}

sub validate_options {
    my ( $gene, $list, $range, $gff_file, $out_dir, $help ) = @_;

    my @errors;

    my @filter_methods = grep { defined $_ } ( $gene, $list, $range );
    push @errors, "Specify ONE filter method: '--gene', '--list', '--range'"
        unless scalar @filter_methods == 1;

    if ( defined $range ) {
        push @errors,
            "When using '--range', must specify a GFF3 annotation file"
            unless defined $gff_file;
    }

    my $provean_out_dir = "$out_dir/pro";
    push @errors, "PROVEAN output directory '$provean_out_dir' not found"
        unless -e $provean_out_dir && -d $provean_out_dir;

    do_or_die( $help, \@errors );
}
