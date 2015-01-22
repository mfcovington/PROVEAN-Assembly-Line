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

use FindBin;
use lib "$FindBin::Bin";
use provean_assembly_line;

my ( $gene, $gene_list_file, $range, $gff_file, $help, $verbose );
my $out_dir = '.';
my $options = GetOptions(
    "gene=s"     => \$gene,
    "list=s"     => \$gene_list_file,
    "range=s"    => \$range,
    "gff_file=s" => \$gff_file,
    "out_dir=s"  => \$out_dir,
    "help"       => \$help,
    "verbose"    => \$verbose,
);

validate_options( $gene, $gene_list_file, $range, $gff_file, $out_dir, $help );

my $gene_list = {};
if ( defined $gene ) {
    $$gene_list{$gene}++;
}
elsif ( defined $gene_list_file ) {
    get_genes_in_list( $gene_list, $gene_list_file );
}
elsif ( defined $range ) {
    get_genes_in_range( $gene_list, $range, $gff_file );
}

my $subs = {};
get_nonsense_subs( 'early', $subs, $gene_list, $out_dir );
get_nonsense_subs( 'late',  $subs, $gene_list, $out_dir );
get_missense_subs( $subs, $gene_list, $out_dir );

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
    my ( $gene_list, $gene_list_file ) = shift;

    open my $gene_list_fh, "<", $gene_list_file;
    for (<$gene_list_fh>) {
        chomp;
        $$gene_list{$_}++;
    }
    close $gene_list_fh;
}

sub get_genes_in_range {
    my ( $gene_list, $range, $gff_file ) = @_;
    my ( $chr, $start, $end ) = split /[:-]/, $range;

    my $genes = get_gene_models($gff_file);

    for my $mrna ( sort keys %{ $$genes{$chr} } ) {
        my $mrna_start = $$genes{$chr}{$mrna}{cds}->[0]->{start};
        my $mrna_end   = $$genes{$chr}{$mrna}{cds}->[-1]->{end};
        $$gene_list{$mrna}++ if $mrna_end >= $start && $mrna_start <= $end;
    }
}

sub get_missense_subs {
    my ( $subs, $gene_list, $out_dir ) = @_;

    for my $gene ( sort keys %$gene_list ) {
        my $provean_file = "$out_dir/pro/$gene.pro";
        my ( $cluster_count, $sequence_count, $provean_scores )
            = parse_provean($provean_file);

        if ( scalar @$provean_scores ) {
            $$subs{$gene}{subs}{clusters} = $cluster_count;
            $$subs{$gene}{subs}{sequences} = $sequence_count;
            $$subs{$gene}{subs}{scores} = $provean_scores;
        }
    }
}

sub get_nonsense_subs {
    my ( $early_or_late, $subs, $gene_list, $out_dir ) = @_;

    open my $stop_fh, "<", "$out_dir/stop/$early_or_late.stop";
    for (<$stop_fh>) {
        chomp;
        my ( $seq_id, @stops ) = split /[\t,]/;
        next unless exists $$gene_list{$seq_id};
        $$subs{$seq_id}{$early_or_late} = \@stops;
    }
    close $stop_fh;
}

sub parse_provean {
    my $provean_file = shift;
    my $cluster_count;
    my $sequence_count;
    my @provean_scores;

    if ( -e $provean_file ) {
        open my $provean_fh, "<", $provean_file;

        my $provean_version = <$provean_fh>;
        validate_provean_version($provean_version);

        for (<$provean_fh>) {
            if (/^\[/) {    # Skip timestamps
                next
            }
            elsif (/# Number of clusters:\t(\d+)/) {
                $cluster_count = $1;
            }
            elsif (/# Number of supporting sequences used:\t(\d+)/) {
                $sequence_count = $1;
            }
            elsif (/^(.+)\t(-?[\d.]+)$/) {
                my $aa_sub = $1;
                my $score  = $2;
                push @provean_scores, [ $aa_sub, $score ];
            }
        }
        close $provean_fh;
    }

    return $cluster_count, $sequence_count, \@provean_scores;
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

sub validate_provean_version {
    my $provean_version = shift;
    my $expected_provean_version = "PROVEAN v1.1 output";
    warn
        "Warning for $gene: Expected $expected_provean_version but found $provean_version"
        if $provean_version !~ /^## $expected_provean_version ##$/;
}

sub write_filtered_results {
    my ( $filtered_out_file ) = @_;

    open my $filtered_out_fh, ">", $filtered_out_file;
    say $filtered_out_fh join "\t", 'gene', 'missense', 'nonsense', 'late stop';
    close $filtered_out_fh;
}
