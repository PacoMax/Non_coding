#!/usr/bin/perl

#=head1 non_cod.pl v.2.0

#=head1 Created

#<2019 - fgonzale>

#=head1 The program gets count table with variables related to coding and noncoding DNA from gfbb files of a directory.

#	The user has to write:

#		./non_cod.pl
		
#		inside the directory

#=cut

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils 'pairwise';
use strict;
use Bio::SeqIO;
use Bio::Location::SplitLocationI;
#my %opts = ();


#GetOptions (\%opts,
#		   'h|help'); 
#if(($opts{'h'}) || (scalar(keys(%opts)) == 0)){
 
#  &PrintHelp();
#}

my @files = glob ("*gbff_*.gbff");

for (my $i = 0; $i <= $#files; $i++){
if (-z $files[$i]){
}
else{
my $ID=$files[$i];
print ("FTP\tGenome_size\tDNA_Cod\tDNA_func\tNum_CDS\tCDS_bp\tNum_tRNA\ttRNA_bp\tNum_rRNA\trRNA_bp\tNum_ncRNA\tncRNA_bp\n");
print ("$ID\t");
open MIFICH, '<', $files[$i];

my $line= <MIFICH>;
if ($line=~/.*\s(\d+)\sbp\s*/){
	print ("$1\t");
	}
close MIFICH;
extraegb($files[$i],my $size=$1);
}
}

sub extraegb {
my $fileGB = $_[0];
my $size= $_[1];
my $Gbfile    = 'not found';
my $locus_tag = 'not found';
my $proteinf  = 'not found';
my $n = 0;
my $f = 0;
my @start_seq = ();
my @end_seq = ();
my @start_seq_fun = ();
my @end_seq_fun = ();
my $n_CDS=0;
my $n_ncRNA=0;
my $n_tRNA=0;
my $n_rRNA=0;
my @end_ncRNA=();
my @start_ncRNA=();
my @end_tRNA=();
my @start_tRNA=();
my @end_rRNA=();
my @start_rRNA=();
my $seqio_object = Bio::SeqIO->new(-file => $fileGB);


	while (my $seq_object = $seqio_object->next_seq){

		for my $feat_object ($seq_object->get_SeqFeatures) {  
			for my $location ($feat_object->location) {
				my $what_is = $feat_object->primary_tag;
				my $start=$location->start;
				my $end=$location->end;
				if($what_is eq 'ncRNA'){
					if ( $location->isa('Bio::Location::SplitLocationI'))  {
						for my $subloc ( $location->sub_Location ) {
							$start_ncRNA[$n_ncRNA]=$subloc->start; 
							$end_ncRNA[$n_ncRNA]=$subloc->end;
							$start_seq_fun[$f]=$start;
							$end_seq_fun[$f]=$end;
						}
					}
					else{
						$start_ncRNA[$n_ncRNA]=$start;
						$end_ncRNA[$n_ncRNA]=$end;
						$start_seq_fun[$f]=$start;
						$end_seq_fun[$f]=$end;
					}
					$n_ncRNA+=1;
					$f+=1;
					}
				if($what_is eq 'rRNA'){
					if ( $location->isa('Bio::Location::SplitLocationI'))  {
						for my $subloc ( $location->sub_Location ) {
							$start_rRNA[$n_rRNA]=$subloc->start; 
							$end_rRNA[$n_rRNA]=$subloc->end;
							$start_seq_fun[$f]=$start;
							$end_seq_fun[$f]=$end;
						}
					}
					else{
						$start_rRNA[$n_rRNA]=$start;
						$end_rRNA[$n_rRNA]=$end;
						$start_seq_fun[$f]=$start;
						$end_seq_fun[$f]=$end;
					}
					$n_rRNA+=1;
					$f+=1;   
				}
				if($what_is eq 'tRNA'){
					if ( $location->isa('Bio::Location::SplitLocationI'))  {
						for my $subloc ( $location->sub_Location ) {
							$start_tRNA[$n_tRNA]=$subloc->start; 
							$end_tRNA[$n_tRNA]=$subloc->end;
							$start_seq_fun[$f]=$start;
							$end_seq_fun[$f]=$end;	
						}
					}
					else{
						$start_tRNA[$n_tRNA]=$start;
						$end_tRNA[$n_tRNA]=$end;
						$start_seq_fun[$f]=$start;
						$end_seq_fun[$f]=$end;
					}
					$n_tRNA+=1;
					$f+=1;   
				}
				if($what_is eq 'CDS'){
					$n_CDS+=1;
					if ( $location->isa('Bio::Location::SplitLocationI'))  {
						for my $subloc ( $location->sub_Location ) {
							$start_seq[$n]=$subloc->start;
							$end_seq[$n]=$subloc->end;
							$start_seq_fun[$f]=$start;
							$end_seq_fun[$f]=$end;	
						}
					}
					else{
						$start_seq[$n]=$start;
						$end_seq[$n]=$end;
						$start_seq_fun[$f]=$start;
						$end_seq_fun[$f]=$end;
					}
					$n+=1;
					$f+=1;	
				}
			}
		}
	}
	
	my @idx = sort { $start_seq[$a] <=> $start_seq[$b] } 0 .. $#start_seq;

	@start_seq = @start_seq[@idx];
	@end_seq = @end_seq[@idx];
	
	my @idx_2 = sort { $start_seq_fun[$a] <=> $start_seq_fun[$b] } 0 .. $#start_seq;

	@start_seq_fun = @start_seq_fun[@idx_2];
	@end_seq_fun = @end_seq_fun[@idx_2];
##########

my $cody=0;
my $l=0;
my $k=-1;
my @end_seq_2=();
my @start_seq_2=();
foreach(@start_seq){
	if ($l>0 && $_<=$end_seq[$l-1]){
		if ($end_seq[$l]>$end_seq[$l-1]){
			$end_seq_2[$k]=$end_seq[$l]
		}
	}
	else{
		$k+=1;
		$start_seq_2[$k]=$_;
		$end_seq_2[$k]=$end_seq[$l];
	}
	$l+=1
}


my @size_cody= pairwise{$a-$b} @end_seq_2,@start_seq_2;
	foreach(@size_cody){
		$cody+=$_+1;
	}
	
my $func=0;
$l=0;
$k=-1;
my @end_seq_fun_2=();
my @start_seq_fun_2=();
foreach(@start_seq_fun){
	if ($l>0 && $_<=$end_seq_fun[$l-1]){
		if ($end_seq_fun[$l]>$end_seq_fun[$l-1]){
			$end_seq_fun_2[$k]=$end_seq_fun[$l]
		}
	}
	else{
		$k+=1;
		$start_seq_fun_2[$k]=$_;
		$end_seq_fun_2[$k]=$end_seq_fun[$l];
	}
	$l+=1
}


my @size_fun= pairwise{$a-$b} @end_seq_fun_2,@start_seq_fun_2;
	foreach(@size_fun){
		$func+=$_+1;
	}

###############
	
							
	my @size_cd= pairwise{$a-$b} @end_seq,@start_seq;
	my $cod=0;
	foreach(@size_cd){
		$cod+=$_+1
	}
	my @size_ncRNA= pairwise{$a-$b} @end_ncRNA,@start_ncRNA;
	my $ncRNA=0;
	foreach(@size_ncRNA){
		$ncRNA+=$_+1
	}
	my @size_tRNA= pairwise{$a-$b} @end_tRNA,@start_tRNA;
	my $tRNA=0;
	foreach(@size_tRNA){
		$tRNA+=$_+1
	}
	my @size_rRNA= pairwise{$a-$b} @end_rRNA,@start_rRNA;
	my $rRNA=0;
	foreach(@size_rRNA){
		$rRNA+=$_+1
	}
	my $total_func=$rRNA+$tRNA+$ncRNA+$cod;
	print ("$cody\t$func\t$n_CDS\t$cod\t$n_tRNA\t$tRNA\t$n_rRNA\t$rRNA\t$n_ncRNA\t$ncRNA\n");

}

#sub PrintHelp {
 #  system "pod2text -c $0 ";
  # exit();
#}
