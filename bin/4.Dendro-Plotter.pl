#!/usr/bin/perl

use File::Basename;
my $programName = basename($0);

use Getopt::Long;
my $usage = "
Usage:   $programName -i input_dir -o output_dir [-s scheme_file]

Arguments: -i  Input wgProfiles directory [String]
           -o  Output dendrogram directory [String]
           -s  Selected scheme file [String]  Optional
             default = 'input_dir/scheme/core.scheme'

Example: Dendro-Plotter with core scheme
         $programName -i ./3.wgProfiles -o ./4.DendroPlot\n\n\n";


my $inPath = undef;
my $outPath = undef;
my $userScheme = undef;

die $usage unless GetOptions(
		'i|input_dir=s'          => \$inPath,
        'o|output_dir=s'         => \$outPath,
        's|scheme_file=s'        => \$userScheme)
	&& defined $inPath 
	&& defined $outPath
#	&& defined $userScheme
	&& @ARGV == 0;

$userScheme = "$inPath/scheme/core.scheme" if($userScheme == undef && -e "$inPath/scheme/core.scheme");

open(MATRIX,"<","$inPath/alleleMatrix") || die "\nError in opening profile matrix $inPath\n\n";
close(MATRIX);

mkdir("$outPath", 0755) || die "$!" if(!-e "$outPath");

#exit();


use FindBin;
my $programPath = $FindBin::Bin;

##### copy files #####

system("cp -rf $inPath/contigfiles.list $outPath/");
system("cp -rf $inPath/locusMatrix $outPath/");
system("cp -rf $inPath/alleleMatrix $outPath/");
system("cp -rf $inPath/locus_profile.xls $outPath/");
system("cp -rf $inPath/allele_profile.xls $outPath/");
system("cp -rf $inPath/locusmapping $outPath/") if(-e "$inPath/locusmapping");
system("cp -rf $inPath/panDB_all.fa $outPath/") if(-e "$inPath/panDB_all.fa");
system("cp -rf $userScheme $outPath/") if(length($userScheme) != 0);

##########

my %filesList = ();

my @filesList = `cat $outPath/contigfiles.list`;
my $assemblyNum = @filesList;
foreach my $filesList(@filesList) {
	chomp($filesList);
	my @temp = split(/\t/, $filesList);
	$temp[0] =~ s/\.fa$//g;
	$temp[1] =~ s/\.\w+$//g;
	$filesList{$temp[0]} = $temp[1];
}

##### read userScheme file #####
my %userScheme = ();
my $schemeCount = 0;

if(length($userScheme) != 0) {
	my @scheme = `cat $userScheme`;
	foreach my $scheme(@scheme) {
		next if ($scheme !~ /^SAL/);
		chomp($scheme);
		$userScheme{$scheme} = 1;
		$schemeCount++;
	}
}
##########



##### read allele matrix #####

my @selScheme = ();
my @alleleDB = ();
my @alleleDB2 = ();

my $locusNum = 0;
my $locusCount = 0;
my @alleleMatrix = `cat $outPath/alleleMatrix`;

foreach my $alleleMatrix(@alleleMatrix) {
	chomp($alleleMatrix);
	my @temp = split(/\t/, $alleleMatrix);
	my $locusIdx = $temp[0];
	next if(length($userScheme) != 0 && $userScheme{$locusIdx} != 1);
	$locusCount++;
	$selScheme[$locusCount] = $locusIdx;
	for(my $i=1;$i<=$assemblyNum;$i++) {
#		my $assemblyID = (sprintf "A%08d", $i);
		$alleleDB[$i][$locusCount] = $temp[$i];
		$alleleDB2[$locusCount][$i] = $temp[$i];
	}
}
$locusNum = $locusCount;

##########


##### calculate allele distance matrix #####

my @alleleDist = ();

for(my $i=1;$i<=$assemblyNum;$i++) {
	my $p1 = $alleleDB[$i];
	for(my $j=1;$j<=$assemblyNum;$j++) {
		my $p2 = $alleleDB[$j];
		$alleleDist[$i][$j] = profileDiff($p1, $p2);
#		print("$i\t$j\t$dist\n");
	}
}

#=head1
  ##### print distance matrix .0.matrix #####
my $alleleMatrixHead = "\t";
for(my $i=1;$i<=$assemblyNum;$i++) {
	my $assemblyID = (sprintf "A%08d", $i);
	if($i != $assemblyNum) {
		$alleleMatrixHead = $alleleMatrixHead . $filesList{$assemblyID} . "\t";
	} else {
		$alleleMatrixHead = $alleleMatrixHead . $filesList{$assemblyID} . "\n";
	}
}
my $alleleMatrixOut = $alleleMatrixHead;
for(my $i=1;$i<=$assemblyNum;$i++) {
	my $assemblyID = (sprintf "A%08d", $i);
	$alleleMatrixOut = $alleleMatrixOut . $filesList{$assemblyID} . "\t";
	for(my $j=1;$j<=$assemblyNum;$j++) {
		if($j != $assemblyNum) {
			$alleleMatrixOut = $alleleMatrixOut . $alleleDist[$i][$j] . "\t";
		} else {
			$alleleMatrixOut = $alleleMatrixOut . $alleleDist[$i][$j] . "\n";
		}
	}
}
open(OUT, ">$outPath/alleleDiff.0.matrix");
print OUT ("$alleleMatrixOut");
close(OUT);
  ##########
#=cut

  ##### print distance matrix .0.LDmatx #####
my $alleleMatrixOut = $assemblyNum . "\n";
for(my $i=1;$i<=$assemblyNum;$i++) {
	my $assemblyID = (sprintf "A%08d", $i);
	$alleleMatrixOut = $alleleMatrixOut . $assemblyID . "\t";
	for(my $j=1;$j<=$i;$j++) {
		if($j != $i) {
			$alleleMatrixOut = $alleleMatrixOut . $alleleDist[$i][$j] . "\t";
		} else {
			$alleleMatrixOut = $alleleMatrixOut . $alleleDist[$i][$j] . "\n";
		}
	}
}
open(OUT, ">$outPath/alleleDiff.0.LDmatx");
print OUT ("$alleleMatrixOut");
close(OUT);
  ##########

##########

##### Bootstraping #####
for(my $n=1;$n<=100;$n++) {
  ##### random shuffling #####
	my @alleleDB_tmp = ();

	for(my $i=1;$i<=$locusNum;$i++) {
		my $random = 1 + int rand($locusNum);
		my $z = $alleleDB2[$random];
#		print("$random\t");
		for(my $j=1;$j<=$assemblyNum;$j++) {
			$alleleDB_tmp[$j][$i] = @$z[$j];
#			print("@$z[$j]\t");
		}
#		print("\n");
	}
  ##########

  ##### calculate distance matrix #####
	my @alleleDist_tmp = ();

	for(my $i=1;$i<=$assemblyNum;$i++) {
		my $p1 = $alleleDB_tmp[$i];
		for(my $j=1;$j<=$assemblyNum;$j++) {
			my $p2 = $alleleDB_tmp[$j];
			$alleleDist_tmp[$i][$j] = profileDiff($p1, $p2);
#			print("$i\t$j\t$dist\n");
		}
	}
  ##########


  ##### print distance matrix .x.LDmatx #####
	$alleleMatrixOut = $assemblyNum . "\n";
	for(my $i=1;$i<=$assemblyNum;$i++) {
		my $assemblyID = (sprintf "A%08d", $i);
		$alleleMatrixOut = $alleleMatrixOut . $assemblyID . "\t";
		for(my $j=1;$j<=$i;$j++) {
			if($j != $i) {
				$alleleMatrixOut = $alleleMatrixOut . $alleleDist_tmp[$i][$j] . "\t";
			} else {
				$alleleMatrixOut = $alleleMatrixOut . $alleleDist_tmp[$i][$j] . "\n";
			}
		}
	}
	open(OUT, ">$outPath/alleleDiff.$n.LDmatx");
	print OUT ("$alleleMatrixOut");
	close(OUT);
  ##########
}

  ##### using neighbor to generate UPGMA trees #####
my $myPath = `pwd`;
chomp($myPath);
chdir("$outPath");
for(my $n=0;$n<=100;$n++) {
	system("mv alleleDiff.$n.LDmatx alleleDiff.LDmatx");
	system("/opt/phylip-3.696/exe/neighbor < $programPath/neighbor.cmd > screen.out");
	if($n == 0) {
		system("cp alleleDiff.LDmatx alleleDiff.$n.LDmatx");
		system("cat outtree > main_tree");
		system("rm -rf bootstrapped_trees") if(-e "bootstrapped_trees");
	} else {
		system("cat outtree >> bootstrapped_trees");
	}
	system("rm -rf outtree") if(-e "outtree");
	system("rm -rf outfile") if(-e "outfile");
	system("rm -rf screen.out") if(-e "screen.out");
}
chdir("$myPath");
  ##########

my $dirName = dirname($outPath);

my $schemeName = basename($userScheme);
$schemeName =~ s/\./ /g;

system("perl -I $programPath $programPath/CompareToBootstrap.pl -tree $outPath/main_tree -boot $outPath/bootstrapped_trees > $outPath/new_tree");
system("xvfb-run python $programPath/drawingTree_1.py $outPath '$dirName' '$schemeName \($schemeCount loci\)'");

##########


sub profileDiff {
	my($p1, $p2) = @_;

	my $totDiff = 0;
	my $tot = @$p1;
	for (my $i=1;$i<=$tot;$i++) {
		$totDiff++ if(@$p1[$i] != @$p2[$i]);
	} 
	return $totDiff;
}

