#!/usr/bin/env perl
# Mike Covington
# created: 2015-01-16
#
# Description:
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';

my @aa_sub_file_list = @ARGV;

for my $aa_sub_file (@aa_sub_file_list) {
    ( my $provean_in_file = $aa_sub_file ) =~ s/aa-changes/provean-in/;

    open my $aa_sub_fh, "<", $aa_sub_file;
    open my $provean_in_fh, ">", $provean_in_file;

    <$aa_sub_fh>;
    while (<$aa_sub_fh>) {
        chomp;
        my ( $seq_id, $aa_subs ) = ( split /\t/ )[ 0, 4 ];
        next unless defined $aa_subs;
        for my $aa_sub ( split /,/, $aa_subs ) {
            my ( $ref, $pos, $alt ) = $aa_sub =~ m/([^\d]+)(\d+)([^\d]+)/;
            say $provean_in_fh join ",", $seq_id, $pos, $ref, $alt;
        }
    }

    close $aa_sub_fh;
    close $provean_in_fh;
}
