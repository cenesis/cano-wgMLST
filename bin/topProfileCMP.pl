#!/usr/bin/perl

if (@ARGV < 3) {
	print("\nUsage: topProfileCMP.pl in_path out_path topScheme_file\n\n");
	print("  perl topProfileCMP.pl ./profileCMP ./topSchemeSel ./topSchemeSel/top5.scheme\n\n");
	exit;
}

use FindBin;
my $programPath = $FindBin::Bin;

my $inPath = $ARGV[0];
my $outPath = $ARGV[1];
my $userScheme = $ARGV[2];

##### copy files #####

system("cp -rf $inPath/contigfiles.list $outPath/");
system("cp -rf $inPath/locusMatrix $outPath/");
system("cp -rf $inPath/alleleMatrix $outPath/");
system("cp -rf $inPath/locus_profile.xls $outPath/");
system("cp -rf $inPath/allele_profile.xls $outPath/");
system("cp -rf $inPath/locusmapping $outPath/") if(-e "$inPath/locusmapping");
system("cp -rf $inPath/tree.newick $outPath/wgMLST_tree.newick") if(-e "$inPath/tree.newick");
system("cp -rf $inPath/panDB_all.fa $outPath/") if(-e "$inPath/panDB_all.fa");
system("cp -rf $userScheme $outPath/used.scheme");

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

use File::Basename;
my $dirName = dirname($outPath);

my $schemeName = basename($userScheme);
$schemeName =~ s/\./ /g;

system("perl -I $programPath $programPath/CompareToBootstrap.pl -tree $outPath/main_tree -boot $outPath/bootstrapped_trees > $outPath/new_tree");
system("xvfb-run python $programPath/drawingTree_2.py $outPath '$dirName' '$schemeName \($schemeCount loci\)'");

##########

##### Robinson-Foulds (RF) distance #####

chdir("$inPath");
my $inPath_abs = `pwd`;
chomp($inPath_abs);
chdir("$myPath");

chdir("$outPath");
system("cat $inPath_abs/tree.newick > two_tree");
system("echo '\n' >> two_tree");
system("cat tree.newick >> two_tree");
system("/opt/phylip-3.696/exe/treedist < $programPath/treedist.cmd > screen.out");
chdir("$myPath");

##########


##### draw heatmap ##### # add on 20190314 (draw heatmap)

my @csv = ();
$csv[0][0] = "#Name";

my @list_heat = `cat $outPath/contigfiles.list`;
my @scheme_heat = `cat $outPath/used.scheme`;

my $isolateNum_heat = @list_heat;
my $locusNum_heat = @scheme_heat;


for($i=0;$i<$isolateNum_heat;$i++) {
	chomp($list_heat[$i]);
	my @temp = split(/\t/, $list_heat[$i]);
	$temp[1] =~ s/\.\w+$//g;
	$csv[$i+1][0] = $temp[1];
}

my $locusCount_heat = 0;
foreach my $locus_heat(@scheme_heat) {
	next if($locus_heat !~ /^SAL/);
	$locusCount_heat++;
	chomp($locus_heat);
	my $alleleMatrix_heat = `grep "^$locus_heat" $outPath/alleleMatrix`;
	chomp($alleleMatrix_heat);
	my @temp1 = split(/\t/, $alleleMatrix_heat);
	for($i=1;$i<=$isolateNum_heat;$i++) {
		$csv[$i][$locusCount_heat] = $temp1[$i];
	}

	my $mapping_heat = `grep "^$locus_heat" $outPath/locusmapping`;
	chomp($mapping_heat);
	my @temp2 = split(/\t/, $mapping_heat);
	my $gene_name = $temp2[1];
	$csv[0][$locusCount_heat] = $gene_name;
}

open(OUT, ">$outPath/profiles.csv");
for($i=0;$i<=$isolateNum_heat;$i++) {
	for($j=0;$j<=$locusNum_heat;$j++) {
		if($j == $locusNum_heat) {
			print OUT ("$csv[$i][$j]\n");
		}else {
			print OUT ("$csv[$i][$j]\t");
		}
	}
}
close(OUT);

system("xvfb-run -a python $programPath/drawHeatmap.py $outPath");

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

