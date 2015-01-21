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
my $out_dir = '.';
my $options = GetOptions(
    "gene=s"    => \$gene,
    "list=s"    => \$list,
    "range=s"   => \$range,
    "out_dir=s" => \$out_dir,
    "help"      => \$help,
    "verbose"   => \$verbose,
);

validate_options( $gene, $list, $range, $out_dir, $help );

my $gene_list;
if ( defined $gene ) {
    push @$gene_list, $gene;
}
elsif ( defined $list ) {

}
elsif ( defined $range ) {

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

sub usage {
    return <<EOF;

Usage: $0 [options] <--filter_method> <filter_argument>

Filter Methods (Specify only one):
  -g, --gene       Specify a single gene (--gene Solyc08g065870.2.1)
  -l, --list       Specify a file with a list of genes (one gene per line)
  -r, --range      Specify a genomic range (--range SL2.40ch01:619126-1047293)

Options:
  -o, --out_dir    Output directory [.]
  -v, --verbose    Print filter summary to STDOUT
  -h, --help       Display this usage information

EOF
}

sub validate_options {
    my ( $gene, $list, $range, $out_dir, $help ) = @_;

    my @errors;

    my @filter_methods = grep { defined $_ } ( $gene, $list, $range );
    push @errors, "Specify ONE filter method: '--gene', '--list', '--range'"
        unless scalar @filter_methods == 1;

    my $provean_out_dir = "$out_dir/pro";
    push @errors, "PROVEAN output directory '$provean_out_dir' not found"
        unless -e $provean_out_dir && -d $provean_out_dir;

    do_or_die( $help, \@errors );
}
