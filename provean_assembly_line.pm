use strict;
use warnings;
use autodie;

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

1;
