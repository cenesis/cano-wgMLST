#!/usr/bin/perl

use File::Basename;
my $programName = basename($0);

use Getopt::Long;
my $usage = "
Usage:   $programName -i input_dir -o output_dir [-s scheme_file] [-t top_threshold]

Arguments: -i  Input dendrogram directory [String]
           -o  Output extracted loci directory [String]
           -s  Selected scheme file [String]  Optional
             default = 'input_dir/core.scheme'
           -t  Top threshold [Integer]  Optional
             default = 5

Example: Loci-Extractor
         $programName -i ./4.DendroPlot -o ./5.ExtractedLoci\n\n\n";


my $inPath = undef;
my $outPath = undef;
my $userScheme = undef;
my $topNum = undef;

die $usage unless GetOptions(
		'i|input_dir=s'          => \$inPath,
        'o|output_dir=s'         => \$outPath,
        's|scheme_file=s'        => \$userScheme,
        't|top=i'                => \$topNum)
	&& defined $inPath 
	&& defined $outPath
#	&& defined $userScheme
#	&& defined $topNum
	&& @ARGV == 0;

$userScheme = "$inPath/core.scheme" if($userScheme == undef && -e "$inPath/core.scheme");
$topNum = 5 if($topNum == undef);

open(TREE,"<","$inPath/new_tree") || die "\nError in opening original tree $inPath\n\n";
close(TREE);

mkdir("$outPath", 0755) || die "$!" if(!-e "$outPath");

#exit();


use FindBin;
my $programPath = $FindBin::Bin;

system("python $programPath/treeTraveling.py $inPath/new_tree $userScheme $inPath/alleleMatrix $outPath");

system("python $programPath/featureImportance.py $outPath $topNum");

system("perl $programPath/topProfileCMP.pl $inPath $outPath $outPath/top$topNum.scheme");

