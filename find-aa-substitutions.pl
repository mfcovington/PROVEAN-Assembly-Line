#!/usr/bin/env perl
# Mike Covington
# created: 2014-04-24
#
# Description:
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use Capture::Tiny 'capture_stderr';
use File::Path 'make_path';
use Getopt::Long;
use Parallel::ForkManager;

use FindBin;
use lib "$FindBin::Bin";
use amino_acid_translation;

my ( $coverage, $par1_bam_file, $par2_bam_file, $ref_vs_alt );

my $gff_file = "ITAG2.3_gene_models.gff3";
my $fa_file  = "S_lycopersicum_chromosomes.2.40.fa";
my $par1_id  = 'M82_n05';
my $par2_id  = 'PEN';
my $threads  = 3;
my $out_dir  = ".";

my $options = GetOptions(
    "gff_file=s"      => \$gff_file,
    "fa_file=s"       => \$fa_file,
    "par1_id=s"       => \$par1_id,
    "par2_id=s"       => \$par2_id,
    "par1_bam_file=s" => \$par1_bam_file,
    "par2_bam_file=s" => \$par2_bam_file,
    "threads=i"       => \$threads,
    "out_dir=s"       => \$out_dir,
    "coverage"        => \$coverage,
    "ref_vs_alt"      => \$ref_vs_alt,
);

my $snp_file_list = [@ARGV];

validate_options( $coverage, $par1_bam_file, $par2_bam_file, $snp_file_list );

my $genes = get_gene_models($gff_file);
my $snps = get_snps( $snp_file_list, $par1_id, $par2_id, $ref_vs_alt );

make_path $out_dir;

my $pm = new Parallel::ForkManager($threads);
for my $seqid ( sort keys %$genes ) {
    next unless exists $$snps{$seqid};
    my $pid = $pm->start and next;

    open my $aa_change_fh, ">", "$out_dir/aa-changes.$seqid";
    write_header( $aa_change_fh, $coverage );

    for my $mrna ( sort keys %{ $$genes{$seqid} } ) {
        my $mrna_start = $$genes{$seqid}{$mrna}{cds}->[0]->{start};
        my $mrna_end   = $$genes{$seqid}{$mrna}{cds}->[-1]->{end};

        filter_snps_by_gene( $snps, $genes, $seqid, $mrna, $mrna_start,
            $mrna_end );

        my ( $par1_seq, $par2_seq )
            = get_seq( $fa_file, $seqid, $mrna_start, $mrna_end,
            $$genes{$seqid}{$mrna}{snps},
            $$snps{$seqid} );

        my ( $par1_cds, $par2_cds, $total_cds_length, $total_cds_coverage )
            = splice_exons( $par1_seq, $par2_seq, $genes, $seqid, $mrna,
            $mrna_start, $coverage );

        if ( $$genes{$seqid}{$mrna}{strand} eq '-' ) {
            reverse_complement( \$par1_cds );
            reverse_complement( \$par2_cds );
        }

        my $par1_protein = translate($par1_cds);
        my $par2_protein = translate($par2_cds);

        my $aa_changes = get_aa_changes( $par1_protein, $par2_protein );

        my $snp_count       = scalar @{ $$genes{$seqid}{$mrna}{snps} };
        my $aa_change_count = scalar @$aa_changes;
        my $change_summary  = join ",", @$aa_changes;

        write_result( $aa_change_fh, $total_cds_coverage, $mrna,
            $total_cds_length, $snp_count, $aa_change_count,
            $change_summary );
    }
    close $aa_change_fh;
    $pm->finish;
}
$pm->wait_all_children;

exit;

sub validate_options {
    my ( $coverage, $par1_bam_file, $par2_bam_file, $snp_file_list ) = @_;

    if ($coverage) {
        die
            "ERROR: --par1_bam_file and --par2_bam_file must be specified when using --coverage\n"
            unless defined $par1_bam_file && defined $par2_bam_file;
    }

    die "ERROR: No SNP file(s) specified\n" if scalar @$snp_file_list == 0;
}

sub get_gene_models {
    my $gff_file = shift;

    my %genes;
    open my $gff_fh, "<", $gff_file;
    while (<$gff_fh>) {
        next if /^#/;
        chomp;

        # GFF Columns:
        # 0:seqid 1:source 2:type 3:start 4:end
        # 5:score 6:strand 7:phase 8:attribute
        my ( $seqid, $type, $start, $end, $strand, $phase, $attribute )
            = ( split /\t/ )[ 0, 2 .. 4, 6 .. 8 ];

        next unless $type =~ /CDS/;

        my ($mrna) = $attribute =~ /Parent=mRNA:([^;]+);/;
        $genes{$seqid}{$mrna}{strand} = $strand;
        push @{ $genes{$seqid}{$mrna}{cds} },
            { start => $start, end => $end, phase => $phase };
    }
    close $gff_fh;

    return \%genes;
}

sub get_snps {
    my ( $snp_file_list, $par1_id, $par2_id, $ref_vs_alt ) = @_;

    my %snps;
    for my $snp_file (@$snp_file_list) {
        open my $snp_fh, "<", $snp_file;
        <$snp_fh>;
        while (<$snp_fh>) {
            my ( $seqid, $pos, $ref, $alt, $alt_parent ) = split /\t/;
            next if $ref =~ /INS/;
            next if $alt =~ /del/;
            if ($ref_vs_alt) {
                $snps{$seqid}{$pos}{par1} = $ref;
                $snps{$seqid}{$pos}{par2} = $alt;
            }
            else {
                $snps{$seqid}{$pos}{par1}
                    = $alt_parent eq $par1_id ? $alt : $ref;
                $snps{$seqid}{$pos}{par2}
                    = $alt_parent eq $par2_id ? $alt : $ref;
            }
        }
        close $snp_fh;
    }

    return \%snps;

    # SNPs data structure:
    # \ {
    #     SEQUENCE ID   {
    #         POSITION   {
    #             par1   PARENT1 ALLELE,
    #             par2   PARENT2 ALLELE
    #         },
    #         POSITION   {
    #             par1   PARENT1 ALLELE,
    #             par2   PARENT2 ALLELE
    #         },
    #         ...
    #     },
    #     ...
    # }
    #
    # Example of hash reference that is returned:
    # \ {
    #     SL2.40ch01   {
    #         13020   {
    #             par1   "C",
    #             par2   "T"
    #         },
    #         13042   {
    #             par1   "C",
    #             par2   "T"
    #         },
    #         ...
    #     },
    #     SL2.40ch02   {
    #         583183    {
    #             par1   "T",
    #             par2   "G"
    #         },
    #         583335    {
    #             par1   "T",
    #             par2   "C"
    #         },
    #         ...
    #     },
    #     ...
    # }
}

sub write_header {
    my ( $aa_change_fh, $coverage ) = @_;

    my @header;
    if ($coverage) {
        @header = (
            'gene',          'length',
            'par1_coverage', 'par2_coverage',
            'snp_count',     'aa_substitution_count',
            'aa_substitutions'
        );
    }
    else {
        @header = (
            'gene',      'length',
            'snp_count', 'aa_substitution_count',
            'aa_substitutions'
        );
    }
    say $aa_change_fh join "\t", @header;
}

sub write_result {
    my ( $aa_change_fh, $total_cds_coverage, $mrna, $total_cds_length,
        $snp_count, $aa_change_count, $change_summary )
        = @_;

    my @results;
    if ($coverage) {
        my $par1_coverage = sprintf "%.1f",
            $$total_cds_coverage{par1} / $total_cds_length;
        my $par2_coverage = sprintf "%.1f",
            $$total_cds_coverage{par2} / $total_cds_length;
        @results = (
            $mrna, $total_cds_length, $par1_coverage, $par2_coverage,
            $snp_count, $aa_change_count, $change_summary
        );
    }
    else {
        @results = (
            $mrna, $total_cds_length, $snp_count, $aa_change_count,
            $change_summary
        );
    }
    say $aa_change_fh join "\t", @results;
}

sub filter_snps_by_gene {
    my ( $snps, $genes, $seqid, $mrna, $mrna_start, $mrna_end ) = @_;

    @{ $$genes{$seqid}{$mrna}{snps} } = ();
    for my $pos ( sort { $a <=> $b } keys %{ $$snps{$seqid} } ) {
        push @{ $$genes{$seqid}{$mrna}{snps} }, $pos
            if ( $pos >= $mrna_start && $pos <= $mrna_end );
    }
}

sub get_seq {
    my ( $fa_file, $seqid, $start, $end, $mrna_snps, $chr_snps ) = @_;

    my ( $header, @seq ) = `samtools faidx $fa_file $seqid:$start-$end`;
    chomp @seq;

    my $par1_seq = join "", @seq;
    my $par2_seq = $par1_seq;

    for my $pos ( @{$mrna_snps} ) {
        my $par1_allele = $$chr_snps{$pos}{par1};
        my $par2_allele = $$chr_snps{$pos}{par2};
        my $offset  = $pos - $start;
        substr $par1_seq, $offset, 1, $par1_allele;
        substr $par2_seq, $offset, 1, $par2_allele;
    }

    return $par1_seq, $par2_seq;
}

sub splice_exons {
    my ( $par1_seq, $par2_seq, $genes, $seqid, $mrna, $mrna_start, $coverage )
        = @_;

    my $par1_cds = '';
    my $par2_cds = '';

    my $total_cds_coverage;
    %$total_cds_coverage = ( par1 => 0, par2 => 0 ) if $coverage;

    for my $cds ( @{ $$genes{$seqid}{$mrna}{cds} } ) {
        my $cds_start = $cds->{start};
        my $cds_end   = $cds->{end};

        if ($coverage) {
            my ( $par1_cds_cov, $par2_cds_cov )
                = get_coverage( $seqid, $cds_start, $cds_end, $fa_file,
                $par1_bam_file, $par2_bam_file );
            $$total_cds_coverage{par1} += $par1_cds_cov;
            $$total_cds_coverage{par2} += $par2_cds_cov;
        }

        my $offset = $cds_start - $mrna_start;
        my $length = $cds_end - $cds_start + 1;
        $par1_cds .= substr $par1_seq, $offset, $length;
        $par2_cds .= substr $par2_seq, $offset, $length;
    }

    my $total_cds_length = length $par1_cds;

    return $par1_cds, $par2_cds, $total_cds_length, $total_cds_coverage;
}

sub reverse_complement {    # Assumes no ambiguous codes
    my $seq = shift;
    $$seq = reverse $$seq;
    $$seq =~ tr/ACGTacgt/TGCAtgca/;
}

sub get_aa_changes {
    my ( $par1_protein, $par2_protein ) = @_;

    my @aa_changes = ();
    for my $idx ( 0 .. length($par1_protein) - 1 ) {
        my $par1_aa = substr $par1_protein, $idx, 1;
        my $par2_aa = substr $par2_protein, $idx, 1;
        next if $par1_aa eq $par2_aa;

        my $pos = $idx + 1;
        push @aa_changes, "$par1_aa$pos$par2_aa";
    }

    return \@aa_changes;
}

sub get_coverage {
    my ( $seqid, $start, $end, $fa_file, $par1_bam_file, $par2_bam_file )
        = @_;
    my $mpileup_cmd
        = "samtools mpileup -A -r $seqid:$start-$end -f $fa_file $par1_bam_file $par2_bam_file";

    my $mpileup_fh;
    capture_stderr {    # suppress mpileup output sent to stderr
        open $mpileup_fh, "-|", $mpileup_cmd;
    };
    my $par1_coverage = 0;
    my $par2_coverage = 0;
    while ( my $mpileup_line = <$mpileup_fh> ) {
        my ($par1_cov_with_gaps, $par1_read_bases,
            $par2_cov_with_gaps, $par2_read_bases
        ) = ( split /\t/, $mpileup_line )[ 3, 4, 6, 7 ];

        # nogap coverage == gap coverage minus # of ref skips
        $par1_coverage += $par1_cov_with_gaps;
        $par2_coverage += $par2_cov_with_gaps;
        $par1_coverage-- for $par1_read_bases =~ m/(?<!\^)[<>]/g;
        $par2_coverage-- for $par2_read_bases =~ m/(?<!\^)[<>]/g;
    }
    close $mpileup_fh;
    return $par1_coverage, $par2_coverage;
}
