use strict;
use autodie;
use Getopt::Long;

my ($fexon, $fgene, $out);
GetOptions(
    'e|exon=s' => \$fexon,
    'g|gene=s' => \$fgene,
    'o|out=s' => \$out
);


# my $fexon = "G3-rna_mm10_1_exon.bam";
# my $fgene = "G3-rna_mm10_1_gene.bam";

open EXON, "samtools view $fexon |";
open GENE, "samtools view $fgene |";

my ($total, $assigned, $three_UTR_exon_count, $five_UTR_exon_count, $CDS_exon_count);

print("Calculating alignment distribution...\n");
while((my $exon = <EXON>) && (my $gene = <GENE>)){
    $total++;
    if($gene =~ /Assigned/){
        $assigned++;
        $five_UTR_exon_count++ if($exon =~ /5_UTR_exon/); 
        $three_UTR_exon_count++ if($exon =~ /3_UTR_exon/); 
        $CDS_exon_count++ if($exon =~ /CDS_exon/);
    }
}

close EXON;
close GENE;

my $intergenic_count = $total-$assigned; 
my $intron = $assigned - $five_UTR_exon_count - $three_UTR_exon_count - $CDS_exon_count;

open OUT, ">", $out;
print OUT "CDS_Exons\t$CDS_exon_count\n";
print OUT "5_UTR_Exons\t$five_UTR_exon_count\n";
print OUT "3_UTR_Exons\t$three_UTR_exon_count\n";
print OUT "Introns\t$intron\n";
print OUT "Intergenic\t$intergenic_count\n";
close OUT;


