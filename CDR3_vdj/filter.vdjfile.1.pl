#!/usr/bin/perl -w
#use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use List::Util qw(max min sum maxstr minstr shuffle);
my $programe_dir=basename($0);
my $path=dirname($0);
my $ver    = "1.0";
my $Writer = "hulongfei\@singleronbio.com";
my $Data   = "2020/08/17";
my $BEGIN=time();
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($inputlist,$namelist,$VDJtype,$celltype,$pair,$loc,$patient,$outdir,$prefix,$datatype);
GetOptions(
			"h|?" =>\&help,
			"inputlist:s"=>\$inputlist,
			"namelist:s"=>\$namelist,
			"VDJtype:s"=>\$VDJtype,
			"Datatype:s"=>\$datatype,
			"celltype:s"=>\$celltype,
			"chain_pair:s"=>\$pair,
			"location:s"=>\$loc,
			"patient:s"=>\$patient,
			"outdir:s"=>\$outdir,
			"prefix:s"=>\$prefix,
			) || &help;
&help unless ($inputlist && $namelist && $prefix && $outdir);
`mkdir $outdir` unless (-e $outdir);
$loc||="NO";
$patient||="NO";
$pair||="NO";
$celltype||="NO";
$VDJtype||="TCR";
$datatype||="10X";
print $loc;
sub help
{
	print <<"	Usage End.";
	Description:
	Writer  : $Writer
	Data    : $Data
	Version : $ver
	function: ......
	Usage:
		-inputlist	input filtered_contig_annotations list,split ,;		must be given

		-namelist	namelist, split ,;					must be given

		-VDJtype	vdj type{TCR,BCR}					[TCR]

		-Datatype	data type{10X,sgr}					[10X]

		-celltype	celltype file by scRNA(ACTION!! barcode match VDJ)	[NO]

		-chain_pair	only use tcr or chain pair barcode 			[NO]

		-location	location TCR,BCR marker					[NO]

		-patient	patient  TCR,BCR marker					[NO]

		-outdir		output file dir						must be given

		-prefix		output file prefix					must be given
		
		-h		Help							document
	Usage End.
	exit;
}
print "input:$inputlist\n";
print "namelist:$namelist\n";
my @invdj=split(/,/,$inputlist);
my @name=split(/,/,$namelist);
my @loc;
my @pat;
if($loc ne "NO"){
	print "location information:$loc\n";
	@loc=split(/,/,$loc);	
}
if($patient ne "NO"){
	print "patient information:$patient\n";
	@pat=split(/,/,$patient);
}
unless(@invdj == @name){die;}
my %celltype;
my %pair_filter;
my %barcode_out;
my %barcode_out_file;
my %barcode_celltype;
my %cellnumber;
if( $celltype ne "NO"){
	open(IN,"$celltype") || die $!;
	while(<IN>){
		chomp;
		my @li=split(/\t/,$_);	
		if(@li != 2){
			next;
		}
		$celltype{$li[0]}=$li[1];
	}
	close IN;
}
my $header;
for (my $i=0;$i<@invdj;$i++){
	if($invdj[$i] =~ /gz$/){
		open(IN,"gzip -dc $invdj[$i]|") || die $!;
	}else{
		open(IN,"$invdj[$i]") || die $!;
	}
	$header=<IN>;
	chomp($header);
	if($datatype ne "10X"){
		$header="barcode,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis";
	}
	if($i==0){
		if($loc ne "NO"){
			$header.=",Location";
		}
		if($patient ne "NO"){
			$header.=",Patient";
		}
	}
	while(<IN>){
		chomp;
		my @li=split(/,/,$_);
		if($datatype ne "10X"){
			my @li1=split(/\t/,$_);
			@li=($li1[0],$li1[1],$li1[2],"None",$li1[3],"None","True","True",$li1[4],$li1[5],$li1[6],$li1[6]);
		}
		if($li[-1] eq "None"){next;}
		my $flag=1;
		if($pair ne "NO"){
			if( $VDJtype eq "TCR" ){
				if($datatype eq "10X"){
					if(($li[5] eq "TRA") and (! exists $pair_filter{$li[0]}{"TRA"})){
						$barcode_out{$li[0]}++;
						$pair_filter{$li[0]}{"TRA"}=1;
						$flag=0;
					}
					if(($li[5] eq "TRB") and (! exists $pair_filter{$li[0]}{"TRB"})){
						$barcode_out{$li[0]}++;
						$pair_filter{$li[0]}{"TRB"}=1;
						$flag=0;
					}
				}else{
					if(($li[1] eq "TRA") and (! exists $pair_filter{$li[0]}{"TRA"}) ){
						$barcode_out{$li[0]}++;
						$pair_filter{$li[0]}{"TRA"}=1;
						$flag=0;
					}
					if(($li[1] eq "TRB") and (! exists $pair_filter{$li[0]}{"TRB"})){
						$barcode_out{$li[0]}++;
						$pair_filter{$li[0]}{"TRB"}=1;
						$flag=0;
					}
				}
			}elsif( $VDJtype eq "BCR" ){
				if(($li[5] eq "IGL") and (! exists $pair_filter{$li[0]}{"IGL"})){
					$barcode_out{$li[0]}++;
					$pair_filter{$li[0]}{"IGL"}=1;
					$flag=0;
				}
				if(($li[5] eq "IGH") and (! exists $pair_filter{$li[0]}{"IGH"})){
					$barcode_out{$li[0]}++;
					$pair_filter{$li[0]}{"IGH"}=1;
					$flag=0;
				}
				if(($li[5] eq "IGK") and (! exists $pair_filter{$li[0]}{"IGK"})){
					$barcode_out{$li[0]}++;
					$pair_filter{$li[0]}{"IGK"}=1;
					$flag=0;
				}
			}else{
				print "$VDJtype is error!!!\n";
			}
			if($flag){next;}
		}
		#%barcode_celltype
		if(($pair ne "NO") and  ($barcode_out{$li[0]} >=2)){
			if(exists $celltype{$li[0]}){
				$cellnumber{$celltype{$li[0]}}++;
				$barcode_celltype{$li[0]}=$celltype{$li[0]};
			}
		}
		if($pair eq "NO"){
			if(($celltype ne "NO") and  (exists $celltype{$li[0]})){
				$cellnumber{$celltype{$li[0]}}++;
				$barcode_celltype{$li[0]}=$celltype{$li[0]}
			}
		}
		if($loc ne "NO"){
			push @li,$loc[$i];
		}
		if($patient ne "NO"){
			push @li,$pat[$i];
		}
		if(exists $barcode_out_file{$li[0]}){
			$barcode_out_file{$li[0]}.="\n".join(",",@li);
		}else{
			$barcode_out_file{$li[0]}=join(",",@li);
		}
#		print OUT join(",",@li)."\n";
	}
	close IN;	
}
#open (OUT,"|gzip > $outdir/$prefix.filtered_contig_annotations.csv.gz")|| die $!;
if( $celltype eq "NO"){
	open (OUT,"|gzip > $outdir/$prefix.filtered_contig_annotations.csv.gz")|| die $!;
	print OUT $header."\n";
	if($pair ne "NO"){
		foreach my $bar (keys %barcode_out){
			if( $barcode_out{$bar} >=2 ){
				print OUT $barcode_out_file{$bar}."\n";
			}
		}
	}
	if($pair eq "NO"){
		foreach my $bar (keys %barcode_out_file){
			print OUT $barcode_out_file{$bar}."\n";
		}
	}
	close OUT;
}else{
	foreach my $cell (sort {$cellnumber{$b} <=> $cellnumber{$a}} keys %cellnumber){
		if($cellnumber{$cell} > 10){
			open ($cell,"|gzip > $outdir/$prefix.$cell.filtered_contig_annotations.csv.gz")|| die $!;
			print {$cell} $header."\n";
		}
	}
	foreach my $barout(keys %barcode_out_file){
		if($pair ne "NO"){
			if( $barcode_out{$barout} >=2  and (exists $celltype{$barout}) and ($cellnumber{$celltype{$barout}} > 10) ){
				print {$celltype{$barout}} $barcode_out_file{$barout}."\n";
			}
		}
		if(($pair eq "NO") and  (exists $celltype{$barout}) and ($cellnumber{$celltype{$barout}} > 10)){
			print {$celltype{$barout}} $barcode_out_file{$barout}."\n";
		}
	}
	foreach my $cell (sort {$cellnumber{$b} <=> $cellnumber{$a}} keys %cellnumber){
		if($cellnumber{$cell} > 10){
			close $cell;
		}
	}
}
