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
use Data::Printer;

my $gff_file = "ITAG2.3_gene_models.gff3";
my @snp_file_list = @ARGV;
# my $snp_file = "M82-PEN.polymorphisms/polyDB.SL2.40ch01.nr";
my $fa_file = "S_lycopersicum_chromosomes.2.40.fa";


my %codon_table = (
    AAA => 'K',
    AAC => 'N',
    AAG => 'K',
    AAT => 'N',
    ACA => 'T',
    ACC => 'T',
    ACG => 'T',
    ACT => 'T',
    AGA => 'R',
    AGC => 'S',
    AGG => 'R',
    AGT => 'S',
    ATA => 'I',
    ATC => 'I',
    ATG => 'M',
    ATT => 'I',
    CAA => 'Q',
    CAC => 'H',
    CAG => 'Q',
    CAT => 'H',
    CCA => 'P',
    CCC => 'P',
    CCG => 'P',
    CCT => 'P',
    CGA => 'R',
    CGC => 'R',
    CGG => 'R',
    CGT => 'R',
    CTA => 'L',
    CTC => 'L',
    CTG => 'L',
    CTT => 'L',
    GAA => 'E',
    GAC => 'D',
    GAG => 'E',
    GAT => 'D',
    GCA => 'A',
    GCC => 'A',
    GCG => 'A',
    GCT => 'A',
    GGA => 'G',
    GGC => 'G',
    GGG => 'G',
    GGT => 'G',
    GTA => 'V',
    GTC => 'V',
    GTG => 'V',
    GTT => 'V',
    TAA => '-',
    TAC => 'Y',
    TAG => '-',
    TAT => 'Y',
    TCA => 'S',
    TCC => 'S',
    TCG => 'S',
    TCT => 'S',
    TGA => '-',
    TGC => 'C',
    TGG => 'W',
    TGT => 'C',
    TTA => 'L',
    TTC => 'F',
    TTG => 'L',
    TTT => 'F',
);

my %cds;
open my $gff_fh, "<", $gff_file;
while (<$gff_fh>) {
    next if /^#/;
    chomp; # Columns: seqid source type start end score strand phase attribute
    my ( $seqid, $type, $start, $end, $strand, $phase, $attribute ) = (split /\t/)[0, 2..4, 6..8];
    next unless $type =~ /CDS/;
    my ($mrna) = $attribute =~ /Parent=mRNA:([^;]+);/;
    $cds{$seqid}{$mrna}{strand} = $strand;
    push @{$cds{$seqid}{$mrna}{cds}}, {start => $start, end => $end, phase => $phase };
}
close $gff_fh;


my $par1 = 'M82_n05';
my $par2 = 'PEN';
my %snps;
for my $snp_file (@snp_file_list) {

    open my $snp_fh, "<", $snp_file;
    <$snp_fh>;
    while (<$snp_fh>) {
        my ( $seqid, $pos, $ref, $alt, $alt_parent ) = split /\t/;
        next if $ref =~ /INS/;
        next if $alt =~ /del/;
        $snps{$seqid}{$pos}{par1} = $alt_parent eq $par1 ? $alt : $ref;
        $snps{$seqid}{$pos}{par2} = $alt_parent eq $par2 ? $alt : $ref;
    }
    close $snp_fh;

}


# my %aa_subs;
open my $aa_change_fh, ">", "aa-changes";

for my $seqid ( sort keys %cds ) {
    next unless exists $snps{$seqid};
# for my $seqid ( 'SL2.40ch01' ) {
    for my $mrna ( sort keys $cds{$seqid} ) {
    # for my $mrna ( 'Solyc01g005000.2.1' ) {
        # next unless $mrna =~ /Solyc01g005/;
        my $mrna_start = $cds{$seqid}{$mrna}{cds}->[0]->{start};
        my $mrna_end   = $cds{$seqid}{$mrna}{cds}->[-1]->{end};
        # $cds{$seqid}{$mrna}{start} = $mrna_start;
        # $cds{$seqid}{$mrna}{end}   = $mrna_end;
        for my $pos ( sort { $a <=> $b } keys $snps{$seqid} ) {
            push @{ $cds{$seqid}{$mrna}{snps} }, $pos
                if ($pos >= $mrna_start && $pos <= $mrna_end );
        }




        # my $seqid = 'SL2.40ch01';
        # my $mrna = 'Solyc01g005000.2.1';
        # p $cds{$seqid}{$mrna};

        # my $mrna_start = 13024;
        # my $mrna_end = 14948;
        # my $seq = get_seq($fa_file, $seqid, $mrna_start, $mrna_end, $cds{$seqid}{$mrna}{snps});
        my ( $par1_seq, $par2_seq ) = get_seq($fa_file, $seqid, $mrna_start, $mrna_end, $cds{$seqid}{$mrna}{snps}, $snps{$seqid});
        # my ( $par1_seq, $par2_seq ) = adjusted_seqs( $seq, $cds{$seqid}{$mrna}{snps} );

        my $ref_spliced = '';
        my $alt_spliced = '';
        for my $cds (@{$cds{$seqid}{$mrna}{cds}}) {
            my $cds_start = $cds->{start};
            my $cds_end = $cds->{end};
            my $cds_phase = $cds->{phase};


            $ref_spliced .= substr $par1_seq, $cds_start - $mrna_start, $cds_end - $cds_start + 1;
            $alt_spliced .= substr $par2_seq, $cds_start - $mrna_start, $cds_end - $cds_start + 1;
        }

        my $ref_protein = translate($ref_spliced);
        my $alt_protein = translate($alt_spliced);

        # $aa_subs{$mrna} = ();
        my @aa_changes = ();
        for my $idx (0 .. length($ref_protein) - 1) {
            my $ref_aa = substr $ref_protein, $idx, 1;
            my $alt_aa = substr $alt_protein, $idx, 1;
            next if $ref_aa eq $alt_aa;
            # push @{$aa_subs{$mrna}}, "$ref_aa:$alt_aa";
            push @aa_changes, "$ref_aa:$alt_aa";
        }

        my $change_count = scalar @aa_changes;
        my $change_summary = join ",", @aa_changes;
        say $aa_change_fh join "\t", $mrna, $change_count, $change_summary;

    }
}

# open my $aa_change_fh, ">", "aa-changes";
# for my $mrna (sort keys %aa_subs) {
#     my @aa_changes = @{$aa_subs{$mrna}} if $aa_subs{$mrna};
#     my $change_count = scalar @aa_changes;
#     my $change_summary = join ",", @aa_changes;
#     say $aa_change_fh join "\t", $mrna, $change_count, $change_summary;
# }
close $aa_change_fh;

# p %aa_subs;

# my $seqid = 'SL2.40ch01';
# my $mrna = 'Solyc01g005000.2.1';
# p $cds{$seqid}{$mrna};

# my $mrna_start = 13024;
# my $mrna_end = 14948;
# # my $seq = get_seq($fa_file, $seqid, $mrna_start, $mrna_end, $cds{$seqid}{$mrna}{snps});
# my ( $par1_seq, $par2_seq ) = get_seq($fa_file, $seqid, $mrna_start, $mrna_end, $cds{$seqid}{$mrna}{snps}, $snps{$seqid});
# # my ( $par1_seq, $par2_seq ) = adjusted_seqs( $seq, $cds{$seqid}{$mrna}{snps} );

# my $ref_spliced = '';
# my $alt_spliced = '';
# for my $cds (@{$cds{$seqid}{$mrna}{cds}}) {
#     my $cds_start = $cds->{start};
#     my $cds_end = $cds->{end};
#     my $cds_phase = $cds->{phase};


#     $ref_spliced .= substr $par1_seq, $cds_start - $mrna_start, $cds_end - $cds_start + 1;
#     $alt_spliced .= substr $par2_seq, $cds_start - $mrna_start, $cds_end - $cds_start + 1;
# }

# my $ref_protein = translate($ref_spliced);
# my $alt_protein = translate($alt_spliced);

# my %aa_subs;
# for my $idx (0 .. length($ref_protein) - 1) {
#     my $ref_aa = substr $ref_protein, $idx, 1;
#     my $alt_aa = substr $alt_protein, $idx, 1;
#     next if $ref_aa eq $alt_aa;
#     push @{$aa_subs{$mrna}}, "$ref_aa:$alt_aa";
# }
# p %aa_subs;

# for my $cds (@{$cds{$seqid}{$mrna}{cds}}) {
#     my $cds_start = $cds->{start};
#     my $cds_end = $cds->{end};
#     my $cds_phase = $cds->{phase};
#     say "$cds_start:$cds_end:$cds_phase";
#     my @cds_snps
#         = get_cds_snps( $cds_start, $cds_end, $cds{$seqid}{$mrna}{snps} );
#     next if scalar @cds_snps == 0;
#     for (@cds_snps) {
#         my $relative_snp_pos = $_ - $mrna_start;
#         my $snp_phase = $relative_snp_pos % 3;
#         my $offset = $relative_snp_pos - abs($snp_phase - $cds_phase);
#         say "PHASES: $snp_phase:$cds_phase";
#         say $relative_snp_pos;
#         say substr $par1_seq, $offset, 3;
#         say substr $par2_seq, $offset, 3;

#     }
#     say "SNPS: @cds_snps";
# }



# say $par1_seq;
# say length $par1_seq;
# say $par2_seq;
# say length $par2_seq;

sub get_seq {
    my ( $fa_file, $seqid, $start, $end, $mrna_snps, $chr_snps ) = @_;
    my ( $header, @seq ) = `samtools faidx $fa_file $seqid:$start-$end`;
    chomp @seq;
    my $par1_seq = join "", @seq;
    # say $par1_seq;
    my $par2_seq = $par1_seq;
    for my $pos (@{$mrna_snps}) {
        my $par1 = $$chr_snps{$pos}{par1};
        my $par2 = $$chr_snps{$pos}{par2};
        my $offset = $pos - $start;
        # say "OFF: $offset.$par1.$par2";
        substr $par1_seq, $offset, 1, $par1;
        substr $par2_seq, $offset, 1, $par2;
        # say "$pos.$par2.$par1";
        # exit;
    }
    return $par1_seq, $par2_seq;
}

# sub get_seq {
#     my ( $fa_file, $seqid, $start, $end ) = @_;
#     my ( $header, @seq ) = `samtools faidx $fa_file $seqid:$start-$end`;
#     chomp @seq;
#     return join "", @seq;
# }

sub get_cds_snps {
    my ( $cds_start, $cds_end, $mrna_snps ) = @_;
    return grep { $_ >= $cds_start && $_ <= $cds_end } @$mrna_snps;
}

sub translate {
    my $nt_seq = shift;
    $nt_seq =~ tr/acgt/ACGT/;

    my $pos = 0;
    my $aa_seq;
    while ( $pos < length($nt_seq) - 2 ) {
        my $codon = substr $nt_seq, $pos, 3;
        my $amino_acid = $codon_table{$codon};
        if ( defined $amino_acid ) {
            $aa_seq .= $amino_acid;
        }
        else {
            $aa_seq .= "X";
        }
        $pos += 3;
    }
    return $aa_seq;
}


