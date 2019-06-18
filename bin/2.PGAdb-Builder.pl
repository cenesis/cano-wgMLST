#!/usr/bin/perl

use File::Basename;
my $programName = basename($0);

use Getopt::Long;
my $usage = "
Usage:   $programName -i input_dir -o output_dir [-m min_identity] [-p number_of_threads]

Arguments: -i  Input cintig annotated directory [String]
           -o  Output PGAdb directory [String]
           -m  Minimum percentage identity for blastp [Integer]  Optional
             default = 95
           -p  Number of threads [Integer]  Optional
             default = 1

Example: Generate PGAdb using 12 threads
         $programName -i ./1.contigAnn -o ./2.PGAdb -p 12\n\n\n";


my $inPath = undef;
my $outPath = undef;
my $identMin = undef;
my $threads = undef;

die $usage unless GetOptions(
		'i|input_dir=s'          => \$inPath,
        'o|output_dir=s'         => \$outPath,
        'm|min_identity=i'       => \$identMin,
        'p|threads=i'            => \$threads)
	&& defined $inPath 
	&& defined $outPath
#	&& defined $identMin
#	&& defined $threads
	&& @ARGV == 0;

$identMin = 95 if($identMin == undef);
$threads = 1 if($threads == undef);

opendir(DIR, "$inPath/GFF") || die "\nError in opening annotated directory $inPath\n\n";
my @gff_files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;
closedir(DIR);
my $gff_num = @gff_files;
die "\nPlease input at least five annotated files.\n\n" if($gff_num < 5);

mkdir("$outPath", 0755) || die "$!" if(!-e "$outPath");


system("xvfb-run roary -p $threads -i $identMin -r -f $outPath/roary $inPath/GFF/*.gff");
system("pdftk $outPath/roary/Rplots.pdf burst output $outPath/roary/Rplot_%02d.pdf");
system("gs -dBATCH -dNOPAUSE -sDEVICE=jpeg -r300 -sOutputFile=$outPath/roary/Rplot_03.jpg $outPath/roary/Rplot_03.pdf");
system("rm -rf doc_data.txt");

##### read all FFN sequences #####

my %ffnDB = ();
my %assemblySeqNum = ();

my @ffnList = `ls $inPath/FFN`;
my $assCount = @ffnList;
foreach my $ffnList(@ffnList) {
	chomp($ffnList);
	my $assName = $ffnList;
	$assName =~ s/\.ffn$//g;
	my @ffnSeq = `cat $inPath/FFN/$ffnList`;
	my $seqLabel = "";
	foreach my $ffnSeq(@ffnSeq) {
		chomp($ffnSeq);
		$ffnSeq =~ s/\r//g;
		if($ffnSeq =~ /^>/) {
			$assemblySeqNum{$assName}++;
			my @labelTmp = split(/ /, $ffnSeq);
			$seqLabel = $labelTmp[0];
			$seqLabel = substr($seqLabel, 10, 5);
#			$seqLabel =~ s/^>PROKKA_//g;
			$seqLabel =~ s/^0+//g;
			$ffnDB{$assName}{$seqLabel} = ">$assName"."_"."$seqLabel\n";
		} else {
			$ffnDB{$assName}{$seqLabel} = $ffnDB{$assName}{$seqLabel} . $ffnSeq;
		}
	}
}

##########


##### read fixed_input_files #####

my %fixMapping = ();

my @fixedList = `ls $inPath/GFF`;
foreach my $fixedList(@fixedList) {
	chomp($fixedList);
	$fixedList =~ s/\.gff$//g;
	my @fixGFF = `grep "locus_tag=" $inPath/GFF/$fixedList.gff`;
	foreach my $fixGFF(@fixGFF) {
		chomp($fixGFF);
		$fixGFF =~ s/[\r\n]//g;
		my @temp1 = split(/\t/, $fixGFF);
		my @temp2 = split(/\;/, $temp1[8]);
		my $id = "";
		my $locus_tag = "";
		foreach my $temp2(@temp2) {
			if($temp2 =~ /^ID=/) {
				$id = $temp2;
				$id =~ s/^ID=//g;
			}
			if($temp2 =~ /^locus_tag=/) {
				$locus_tag = $temp2;
				$locus_tag =~ s/^locus_tag=//g;
				$locus_tag = substr($locus_tag, 9, 5);
				$locus_tag =~ s/^0+//g;
			}
		}
		$fixMapping{$id} = $locus_tag;

	}
}

##########


##### read gene_presence_absence.csv file #####

my %csvDB = ();
my $locusCount = 0;
my @csvHeader = ();
my @locusfilesList = ();
my $paralog = "";

my @csv = `cat $outPath/roary/gene_presence_absence.csv`;
foreach my $csv(@csv) {
	chomp($csv);
	$csv =~ s/\r//g;
	$csv =~ s/^\"//g;
	$csv =~ s/\"$//g;
	my @temp = split(/\"\,\"/, $csv);

	if($csv =~ /^Gene/) {
		@csvHeader = split(/\"\,\"/, $csv);
		next;	
	}

	if ($temp[3] ne $temp[4]){
#		print("$csv\n");
		$paralog = $paralog . "$temp[0]\t$temp[3]\t$temp[4]\t$temp[2]\n";
		next;
	}

	$locusCount++;
	my $locusID = (sprintf "SAL%07d", $locusCount);
	push(@locusfilesList, $locusID);

	$csvDB{$locusID}{"gene"} = $temp[0];
	$csvDB{$locusID}{"anno"} = $temp[2];
	$csvDB{$locusID}{"nIso"} = $temp[3];
	$csvDB{$locusID}{"nSeq"} = $temp[4];

	for(my $i=14;$i<$assCount+14;$i++) {
		if($i == 14) {
			my @tmp = split(/_/, $temp[$i]);
			$tmp[1] =~ s/^0+//g;
			$csvDB{$locusID}{$csvHeader[$i]} = $tmp[1];
		} else {
			$csvDB{$locusID}{$csvHeader[$i]} = $fixMapping{$temp[$i]};
		}
	}
}

open(OUT1, ">$outPath/locusfiles.list");
foreach my $locusfilesList(@locusfilesList) {
	print OUT1 ("$locusfilesList\n");
}
close(OUT1);

open(OUT2, ">$outPath/locusmapping");
print OUT2 ("Locus\tGene\tNo. isolates\tNo. sequences\tAnnotation\n");
foreach my $locusfilesList(@locusfilesList) {
	print OUT2 ("$locusfilesList\t$csvDB{$locusfilesList}{'gene'}\t$csvDB{$locusfilesList}{'nIso'}\t$csvDB{$locusfilesList}{'nSeq'}\t$csvDB{$locusfilesList}{'anno'}\n");
}
close(OUT2);

open(OUT3, ">$outPath/paralog.txt");
print OUT3 ("$paralog");
close(OUT3);


##########


##### write locus files #####

my %panRefSeq = ();

mkdir("$outPath/locusfiles", 0755) || die "$!" if(!-e "$outPath/locusfiles");
foreach my $locusfilesList(@locusfilesList) {
	open(OUT, ">$outPath/locusfiles/$locusfilesList.tmp");
	for(my $i=14;$i<$assCount+14;$i++) {
		if($csvDB{$locusfilesList}{$csvHeader[$i]} ne "") {
			print OUT ("$ffnDB{$csvHeader[$i]}{$csvDB{$locusfilesList}{$csvHeader[$i]}}\n");
		}
	}
	close(OUT);
	system("fastx_collapser -i $outPath/locusfiles/$locusfilesList.tmp -o $outPath/locusfiles/$locusfilesList");
	
	##### renumber alleles #####
	my @locusFiles = `cat $outPath/locusfiles/$locusfilesList`;
	my $locusFilesOut = "";
	my $panRefSeqLen = -1;
	my $panRefSeq = "";
	foreach $locusFiles(@locusFiles) {
		if($locusFiles =~ /^>/) {
			chomp($locusFiles);
			my @temp = split(/\-/, $locusFiles);
			$temp[0] =~ s/^>//g;
			$locusFiles = ">" . $locusfilesList . "::" . $temp[0] . "\n";
			$locusFilesOut = $locusFilesOut . $locusFiles;
		} else {
			$locusFilesOut = $locusFilesOut . $locusFiles;
			if(length($locusFiles) > $panRefSeqLen) {
				$panRefSeqLen = length($locusFiles);
				$panRefSeq = ">" . $locusfilesList . "\n" . $locusFiles;
			}
		}
	}
	open(OUT, ">$outPath/locusfiles/$locusfilesList");
	print OUT ("$locusFilesOut");
	close(OUT);
	##########

	$panRefSeq{$locusfilesList} = $panRefSeq;

}
##########

##### make Pan References #####

open(OUT, ">$outPath/panRefSeq.fa");
foreach my $locusfilesList(@locusfilesList) {
	print OUT ("$panRefSeq{$locusfilesList}");
}
close(OUT);

##########


##### write scheme files (scheme.txt, core.txt, pan.txt) #####
mkdir("$outPath/scheme", 0755) || die "$!" if(!-e "$outPath/scheme");

my $core100List = "core100\n";
my $core99List = "core99\n";
my $core95List = "core95\n";
my $core90List = "core90\n";
my $core80List = "core80\n";
my $core70List = "core70\n";
my $core60List = "core60\n";
my $core50List = "core50\n";
my $core30List = "core30\n";
my $core10List = "core10\n";

my $coreList = "core\n";
my $accessoryList = "accessory\n";
my $panList = "pan\n";
my $uniqueList = "unique\n";

foreach my $locusfilesList(@locusfilesList) {
	my $occ = $csvDB{$locusfilesList}{'nIso'} / $assCount;

	if($occ >= 1.0) {
		$core100List = $core100List . $locusfilesList . "\n";
	}
	if($occ >= 0.99) {
		$core99List = $core99List . $locusfilesList . "\n";
	}
	if($occ >= 0.95) {
		$core95List = $core95List . $locusfilesList . "\n";
	}
	if($occ >= 0.90) {
		$core90List = $core90List . $locusfilesList . "\n";
	}
	if($occ >= 0.80) {
		$core80List = $core80List . $locusfilesList . "\n";
	}
	if($occ >= 0.70) {
		$core70List = $core70List . $locusfilesList . "\n";
	}
	if($occ >= 0.60) {
		$core60List = $core60List . $locusfilesList . "\n";
	}
	if($occ >= 0.50) {
		$core50List = $core50List . $locusfilesList . "\n";
	}
	if($occ >= 0.30) {
		$core30List = $core30List . $locusfilesList . "\n";
	}
	if($occ >= 0.10) {
		$core10List = $core10List . $locusfilesList . "\n";
	}

	if($occ >= 0.95 && $occ <= 1) {
		$coreList = $coreList . $locusfilesList . "\n";
	} elsif ($csvDB{$locusfilesList}{'nIso'} != 1) {
		$accessoryList = $accessoryList . $locusfilesList . "\n";
	}

	if($csvDB{$locusfilesList}{'nIso'} == 1) {
		$uniqueList = $uniqueList . $locusfilesList . "\n";
	}
	$panList = $panList . $locusfilesList . "\n";
}

open(OUT, ">$outPath/scheme/core100.scheme");
print OUT ("$core100List");
close(OUT);
open(OUT, ">$outPath/scheme/core99.scheme");
print OUT ("$core99List");
close(OUT);
open(OUT, ">$outPath/scheme/core95.scheme");
print OUT ("$core95List");
close(OUT);
open(OUT, ">$outPath/scheme/core90.scheme");
print OUT ("$core90List");
close(OUT);
open(OUT, ">$outPath/scheme/core80.scheme");
print OUT ("$core80List");
close(OUT);
open(OUT, ">$outPath/scheme/core70.scheme");
print OUT ("$core70List");
close(OUT);
open(OUT, ">$outPath/scheme/core60.scheme");
print OUT ("$core60List");
close(OUT);
open(OUT, ">$outPath/scheme/core50.scheme");
print OUT ("$core50List");
close(OUT);
open(OUT, ">$outPath/scheme/core30.scheme");
print OUT ("$core30List");
close(OUT);
open(OUT, ">$outPath/scheme/core10.scheme");
print OUT ("$core10List");
close(OUT);

open(OUT, ">$outPath/scheme/core.scheme");
print OUT ("$coreList");
close(OUT);

open(OUT, ">$outPath/scheme/accessory.scheme");
print OUT ("$accessoryList");
close(OUT);

open(OUT, ">$outPath/scheme/pan.scheme");
print OUT ("$panList");
close(OUT);

open(OUT, ">$outPath/scheme/unique.scheme");
print OUT ("$uniqueList");
close(OUT);

system("rm -rf $outPath/locusfiles/*.tmp");

##########


