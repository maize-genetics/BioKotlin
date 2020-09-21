use strict;
use warnings;
use List::Util qw[min max];
my ($pfam_pos_file,$gff,$pfam_genome_pos)= @ARGV;
 open OUT,">","$pfam_genome_pos";
open ALL_PEP_PFAM,"$pfam_pos_file";
my @pep_pfam=<ALL_PEP_PFAM>;

print OUT "UniProtIDt\tZmProteinID\tpep_chr\tpep_start\tpep_end\tdomain_id\tdomainScore\tdomainStart\tdomainEnd\tdomain_genome_start\tdomain_genome_end\n";
my %hash_pep;
open GFF,"$gff";
my @gff=<GFF>;
my $all_gffs="@gff";
my @gff_infor=split/###/,$all_gffs;

foreach my $one_gene_gff(@gff_infor) {
    my @one_gene_all_cds;
    my @one_gene=split /\n/,$one_gene_gff;
    my $aa_start =0;
    my $aa_end =0;
    my $pep_id="";
    my $previous_pep="";
    foreach my $one_line(@one_gene){
        chomp($one_gene_gff);
        $one_line =~ s/^ //g;
        if($one_line=~ /ID=CDS:/) {
            my @gene_gff_infor = split /\t/, $one_line;
            my @pep_infor = split /[:;]/, $gene_gff_infor[8];
            $pep_id = $pep_infor[1];
            if ($pep_id eq $previous_pep){
                $aa_start = $aa_end + 1;
            } else {
                $hash_pep{$previous_pep}="@one_gene_all_cds";
                @one_gene_all_cds="";
                $aa_start =0;
            }
            $aa_end = $aa_start + int(($gene_gff_infor[4] - $gene_gff_infor[3] + 1 )/ 3) - 1;
            my $CDS_pep_pos = "$gene_gff_infor[0]\t$gene_gff_infor[2]\t$gene_gff_infor[3]\t$gene_gff_infor[4]\t$aa_start\t$aa_end\t$pep_id\n";
            push @one_gene_all_cds, $CDS_pep_pos;
            $previous_pep=$pep_id;
        }
    }
    $hash_pep{$pep_id}="@one_gene_all_cds";
}
shift(@pep_pfam);
foreach my $one_pep_pfam(@pep_pfam) {
    chomp($one_pep_pfam);
    my @pep_domain_infor = split /\s+/, $one_pep_pfam;
    my $domain_satrt = $pep_domain_infor[4];
    my $domain_end = $pep_domain_infor[5];
    my $domain_id=$pep_domain_infor[2];
    my $gene_pep_id = $pep_domain_infor[1];
    my @CDS_infor =split /\n/, $hash_pep{$gene_pep_id};
    my $domain_genome_start=0;
    my $domain_genome_end=0;
    my $pep_chr;
    my $pep_start=10000000000;
    my $pep_end=0;
    foreach my $CDS (@CDS_infor) {
        $CDS =~ s/ //g;
        my @CDS_pos = split /\s+/, $CDS;
        $pep_chr=$CDS_pos[0];
        if ($CDS_pos[4] <= $domain_satrt && $CDS_pos[5] >= $domain_satrt) {
            $domain_genome_start = $domain_satrt * 3 + $CDS_pos[2];
        }
        if ($CDS_pos[4] <= $domain_end && $CDS_pos[5] >= $domain_end) {
            $domain_genome_end = $domain_end * 3 + $CDS_pos[2];
        }
        $pep_start=min($pep_start,$CDS_pos[2]);
        $pep_end=max($pep_end,$CDS_pos[3]);
    }
     print OUT "$pep_domain_infor[0]\t$gene_pep_id\t$pep_chr\t$pep_start\t$pep_end\t$domain_id\t$domain_satrt\t$domain_end\t$domain_genome_start\t$domain_genome_end\n";
}
