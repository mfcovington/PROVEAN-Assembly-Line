#!/usr/bin/env perl
# Mike Covington
# created: 2015-01-16
#
# Description: Run PROVEAN on results in aa-changes.* files
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use amino_acid_translation;
use File::Path 'make_path';
use Getopt::Long;

my ( $cds_fasta_file, $supporting_set, $verbose );
my $threads = 1;
my $out_dir = '.';

my $options = GetOptions(
    "cds_fasta_file=s" => \$cds_fasta_file,
    "threads=i"        => \$threads,
    "out_dir=s"        => \$out_dir,
    "supporting_set"   => \$supporting_set,
    "verbose"          => \$verbose,
);

my $aa_sub_file_list = [@ARGV];

validate_options( $cds_fasta_file, $threads, $aa_sub_file_list );

make_path "$out_dir/$_" for ( 'fa', 'pro', 'sss', 'var' );

for my $aa_sub_file (@$aa_sub_file_list) {
    open my $aa_sub_fh, "<", $aa_sub_file;
    <$aa_sub_fh>;
    while (<$aa_sub_fh>) {
        chomp;
        my ( $seq_id, $aa_subs ) = ( split /\t/ )[ 0, 4 ];
        next unless defined $aa_subs;

        write_provean_fa_file( $seq_id, $cds_fasta_file, $out_dir );
        write_provean_var_file( $seq_id, $aa_subs, $out_dir );
        run_provean( $seq_id, $out_dir, $threads, $supporting_set, $verbose );
    }
    close $aa_sub_fh;
}

exit;

sub run_provean {
    my ( $seq_id, $out_dir, $threads, $supporting_set, $verbose ) = @_;

    my $provean_cmd = <<EOF;
provean.sh \\
  -q $out_dir/fa/$seq_id.fa \\
  -v $out_dir/var/$seq_id.var \\
  --num_threads $threads \\
EOF

    if ($supporting_set) {
        $provean_cmd .= "  --supporting_set $out_dir/sss/$seq_id.sss";
    }
    else {
        $provean_cmd .= "  --save_supporting_set $out_dir/sss/$seq_id.sss";
    }

    $provean_cmd .= " \\\n  > $out_dir/pro/$seq_id.pro\n";

    say STDERR $seq_id if $verbose;
    system($provean_cmd);
}

sub usage {
    return "<USAGE STATEMENT PLACEHOLDER>\n";
}

sub validate_options {
    my ( $cds_fasta_file, $threads, $aa_sub_file_list ) = @_;

    my @errors;

    if ( defined $cds_fasta_file ) {
        push @errors, "File '$cds_fasta_file' not found"
            if !-e $cds_fasta_file;
    }
    else {
        push @errors, "Must specify '--cds_fasta_file'";
    }

    push @errors, "Option '--threads' must be an integer greater than 0"
        if $threads <= 0;

    push @errors, "Must specify amino acid substitution files"
        if scalar @$aa_sub_file_list == 0;

    if ( scalar @errors > 0 ) {
        my $error_string = join "\n", map {"ERROR: $_"} @errors;
        die join( "\n", usage(), $error_string ), "\n";
    }
}

sub write_provean_fa_file {
    my ( $seq_id, $cds_fasta_file, $out_dir ) = @_;
    my $fa_width = 80;

    open my $cds_fa_in_fh, "-|", "samtools faidx $cds_fasta_file $seq_id";
    my ( $header, @cds_seq ) = <$cds_fa_in_fh>;
    close $cds_fa_in_fh;

    chomp @cds_seq;
    my $aa_seq = translate( join "", @cds_seq );

    open my $provean_fa_out_fh, ">", "$out_dir/fa/$seq_id.fa";
    print $provean_fa_out_fh $header;
    say $provean_fa_out_fh $_ for unpack "(A$fa_width)*", $aa_seq;
    close $provean_fa_out_fh;
}

sub write_provean_var_file {
    my ( $seq_id, $aa_subs, $out_dir ) = @_;

    open my $provean_var_out_fh, ">", "$out_dir/var/$seq_id.var";
    say $provean_var_out_fh $_ for split /,/, $aa_subs;
    close $provean_var_out_fh;
}
