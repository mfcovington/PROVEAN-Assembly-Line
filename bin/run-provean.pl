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
use File::Path 'make_path';
use Getopt::Long;
use Parallel::ForkManager;
use POSIX 'strftime';

use FindBin;
use lib "$FindBin::Bin/../lib";
use amino_acid_translation;

my ( $cds_fasta_file, $supporting_set, $force, $verbose, $help );
my $threads = 1;
my $out_dir = '.';

my $options = GetOptions(
    "cds_fasta_file=s" => \$cds_fasta_file,
    "threads=i"        => \$threads,
    "out_dir=s"        => \$out_dir,
    "supporting_set"   => \$supporting_set,
    "force"            => \$force,
    "verbose"          => \$verbose,
    "help"             => \$help,
);

my $aa_sub_file_list = [@ARGV];

validate_options( $cds_fasta_file, $threads, $out_dir, $aa_sub_file_list,
    $force, $help );

make_path "$out_dir/$_" for ( 'fa', 'pro', 'sss', 'stop', 'var' );

my $pm = new Parallel::ForkManager($threads);
for my $aa_sub_file (@$aa_sub_file_list) {
    open my $aa_sub_fh, "<", $aa_sub_file;
    <$aa_sub_fh>;
    while (<$aa_sub_fh>) {
        chomp;
        my ( $seq_id, $aa_subs ) = ( split /\t/ )[ 0, 4 ];
        next unless defined $aa_subs;

        my $pid = $pm->start and next;
        write_provean_fa_file( $seq_id, $cds_fasta_file, $out_dir );
        write_provean_var_file( $seq_id, $aa_subs, $out_dir );
        run_provean( $seq_id, $out_dir, $supporting_set, $verbose );
        $pm->finish;
    }
    close $aa_sub_fh;
}
$pm->wait_all_children;

exit;

sub run_provean {
    my ( $seq_id, $out_dir, $supporting_set, $verbose ) = @_;

    my $provean_cmd = <<EOF;
provean.sh \\
  -q $out_dir/fa/$seq_id.fa \\
  -v $out_dir/var/$seq_id.var \\
EOF

    if ($supporting_set) {
        $provean_cmd .= "  --supporting_set $out_dir/sss/$seq_id.sss";
    }
    else {
        $provean_cmd .= "  --save_supporting_set $out_dir/sss/$seq_id.sss";
    }

    $provean_cmd .= " \\\n  > $out_dir/pro/$seq_id.pro\n";

    if ($verbose) {
        my $timestamp = strftime "%Y-%m-%d %H:%M:%S", localtime;
        say STDERR "$timestamp - $seq_id" if $verbose;
    }

    if ( -s "$out_dir/var/$seq_id.var" ) {
        system($provean_cmd);
    }
    else {return}    # Don't run if only aa subs were early/late STOP codons
}

sub usage {
    return <<EOF;

Usage: $0 [options] --cds_fasta_file <CDS.fa> <Amino acid change file(s)>

Options:
  -t, --threads           Number of threads to use for PROVEAN [1]
  -o, --out_dir           Output directory [.]
  -s, --supporting_set    Supporting set files already exist
  -v, --verbose           Report current progress of PROVEAN
  -h, --help              Display this usage information

EOF
}

sub validate_options {
    my ( $cds_fasta_file, $threads, $out_dir, $aa_sub_file_list, $force, $help ) = @_;

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

    unless ($force) {
        push @errors,
            "File $out_dir/stop/early.stop already exists (Choose a new --out_dir, remove file, or use --force to overwrite)"
            if -e "$out_dir/stop/early.stop";

        push @errors,
            "File $out_dir/stop/late.stop already exists (Choose a new --out_dir, remove file, or use --force to overwrite)"
            if -e "$out_dir/stop/late.stop";
    }

    push @errors, "PROVEAN not installed in PATH" unless `which provean`;

    if (@errors) {
        my $error_string = join "\n", map {"ERROR: $_"} @errors;
        die usage(), $error_string, "\n\n";
    }
    elsif ($help) {
        die usage();
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
    $aa_seq =~ s/-$//;    # PROVEAN can't deal with '-'

    say STDERR "WARNING: $seq_id contains an unexpected interior STOP codon"
        if $aa_seq =~ /-/;

    open my $provean_fa_out_fh, ">", "$out_dir/fa/$seq_id.fa";
    print $provean_fa_out_fh $header;
    say $provean_fa_out_fh $_ for unpack "(A$fa_width)*", $aa_seq;
    close $provean_fa_out_fh;
}

sub write_provean_var_file {
    my ( $seq_id, $aa_subs, $out_dir ) = @_;
    my @early_stop;

    open my $provean_var_out_fh, ">", "$out_dir/var/$seq_id.var";
    for my $aa_sub ( split /,/, $aa_subs ) {

        if ( $aa_sub =~ /-$/ ) {
            push @early_stop, $aa_sub;
        }
        elsif ( $aa_sub =~ /^-/ ) {
            open my $provean_late_out_fh, ">>", "$out_dir/stop/late.stop";
            say $provean_late_out_fh join "\t", $seq_id, $aa_sub;
            close $provean_late_out_fh;
        }
        else {
            say $provean_var_out_fh $aa_sub;
        }
    }
    close $provean_var_out_fh;

    if (@early_stop) {
        open my $provean_early_out_fh, ">>", "$out_dir/stop/early.stop";
        say $provean_early_out_fh $seq_id, "\t", join ",", @early_stop;
        close $provean_early_out_fh;
    }
}
