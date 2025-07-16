use List::Util qw(max);
use strict;
use autodie;
use Getopt::Long;
use File::Basename; 


# perl ./make_10X.pl -g /home/zeemeeuw/YangLab/ref/hg19/hg19.refGene.gtf -i /home/zeemeeuw/YangLab/Ti-ATAC-seq2/Data/JH-rna/JH-rna_counts.tsv.gz -o /home/zeemeeuw/YangLab/Ti-ATAC-seq2/Data/JH-rna/JH-rna_out.tar.gz
my ($gtf, $infile, $outfile, $min_umi);
$min_umi = 100;
GetOptions(
    'g|gtf=s' => \$gtf,
    'i|in=s' => \$infile,
    'o|out=s' => \$outfile,
    'm|min:i' => \$min_umi
);

BEGIN{$, = "\t";}

# get dirname
my $dir = dirname($outfile);
my $filtfile = $dir."/"."filter.tsv";
my $matrix = $dir."/"."matrix.mtx";
my $barcode = $dir."/"."barcodes.tsv";
my $feature = $dir."/"."features.tsv";
my @genes;
my @cells;
my @counts;
my $total;
my %umi_counts;

if ($infile =~ /\.gz\z/){
    open IN, "zcat $infile | ";
}else{
    open IN, "<", $infile;
}


while(<IN>){
    if($.>1){
        chomp;
        my @F = split /\t/;
        $umi_counts{$F[1]} += $F[2];
    }
}
close IN;


open OUT, ">", $filtfile;
if ($infile =~ /\.gz\z/){
    open IN, "zcat $infile | ";
}else{
    open IN, "<", $infile;
}
while(<IN>){
    if($.>1){
        my @F = split /\t/;
        print OUT $_ if $umi_counts{$F[1]} >= $min_umi;
    }else{print OUT $_;}
}
close IN;
close OUT;

system("gzip -f $filtfile");
my $tmp = $filtfile.".gz";


open IN, "zcat $tmp | ";
while(<IN>){
    if($.>1){
        chomp;
        my @F = split /\t/;
        push(@genes, $F[0]);
        push(@cells, $F[1]);
        push(@counts, $F[2]);
        # $total += $F[2];
        $total += 1;
    }
}
close IN;

# create gene symbols - gene id hash
my %gid;
open GTF, "<", $gtf;
while(<GTF>){
    unless(/^#/){
        chomp;
        my @F = split /\t/;
        if($F[8] =~ /gene_id "(.*?)";.*gene_name "(.*?)";/){$gid{$2} = $1;}
    }
}
close GTF;

my @genes_index = &name_to_index_genes($feature, @genes);
my @cells_index = &name_to_index_cells($barcode, @cells);

open MTX, ">", $matrix;
my @out;
foreach(0..$#counts){push(@out, [$genes_index[$_], $cells_index[$_], $counts[$_]]);}
my $genes_max = max(@genes_index);
my $cells_max = max(@cells_index);
print MTX "%%MatrixMarket\tmatrix\tcoordinate\tinteger\tgeneral\n";
print MTX "%\n";
print MTX "$genes_max\t$cells_max\t$total\n";
my @sorted_out = sort { $b->[0] <=> $a->[0] or $b->[1] <=> $a->[1] or $b->[2] <=> $a->[2] } @out;   
foreach(@sorted_out){print MTX "$_->[0]\t$_->[1]\t$_->[2]\n";}
close MTX;

system("gzip -f $matrix");
system("gzip -f $barcode");
system("gzip -f $feature");
my ($b_barcode, $b_matrix, $b_feature) = (basename($barcode).".gz", basename($matrix).".gz", basename($feature).".gz");
system("tar --remove-files -C $dir -zcvf $outfile $b_matrix $b_barcode $b_feature");

system("rm $tmp");

sub name_to_index_cells{
    my $n = 1;
    my @index;
    my %h;
    my %i;
    my $out = shift @_;
    foreach(@_){
        $h{$_}++;
        unless($h{$_}>1){
            push(@index, $n);
            $i{$_} = $n;
            $n++;
        }else{
            push(@index, $i{$_});
        }
    }
    open OUT, ">", $out;
    foreach(sort { $i{$a} <=> $i{$b} } keys %i){
        print OUT "$_\n";
    }
    close OUT;
    return @index;
}

sub name_to_index_genes{
    my $n = 1;
    my @index;
    my %h;
    my %i;
    my $out = shift @_;
    foreach(@_){
        $h{$_}++;
        unless($h{$_}>1){
            push(@index, $n);
            $i{$_} = $n;
            $n++;
        }else{
            push(@index, $i{$_});
        }
    }
    open OUT, ">", $out;
    foreach(sort { $i{$a} <=> $i{$b} } keys %i){
        print OUT "$gid{$_}\t$_\n";
    }
    close OUT;
    return @index;
}


