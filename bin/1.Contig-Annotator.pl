#!/usr/bin/perl

use File::Basename;
my $programName = basename($0);

use Getopt::Long;
my $usage = "
Usage:    $programName -i input_dir -o output_dir [-p number_of_threads]

Arguments: -i  Input contig directory [String]
           -o  Output annotated contig directory [String]
           -p  Number of threads [Integer]  Optional
             default = 1

Example: Generate annotated contig files using 12 threads
         $programName -i ./contigFiles -o ./1.contigAnn -p 12\n\n\n";


my $inPath = undef;
my $outPath = undef;
my $threads = undef;

die $usage unless GetOptions(
		'i|input_dir=s'          => \$inPath,
        'o|output_dir=s'         => \$outPath,
        'p|threads=i'            => \$threads)
	&& defined $inPath 
	&& defined $outPath
#	&& defined $threads
	&& @ARGV == 0;

$threads = 1 if($threads == undef);

opendir(DIR, "$inPath") || die "\nError in opening contigs directory $inPath\n\n";
my @contig_files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;
closedir(DIR);
my $contig_num = @contig_files;
die "\nPlease input at least five genome contig files.\n\n" if($contig_num < 5);

mkdir("$outPath", 0755) || die "$!" if(!-e "$outPath");


mkdir("$outPath/Assembly", 0755) || die "$!" if(!-e "$outPath/Assembly");
my @list = `ls $inPath`;
my $m = 0;
my $newList = "";
foreach my $list(@list) {
	chomp($list);
	$m++;
	my $num_m = (sprintf "%08d", $m);
	my $n = 0;
	my @temp = `cat $inPath/$list`;
	my $seq = "";
	foreach my $temp(@temp) {
		chomp($temp);
		if($temp =~ /^>/) {
			$n++;
			my $num_n = (sprintf "%08d", $n);
			if($n == 1) {
				$seq = $seq . ">A$num_m\:\:C$num_n" . "\n";
			} else {
				$seq = $seq . "\n>A$num_m\:\:C$num_n" . "\n";
			}
		} else {
			$seq = $seq . $temp;
		}
	}
	$seq = $seq . "\n";

	open(OUT, ">$outPath/Assembly/A$num_m.fa");
	print OUT ("$seq");
	close(OUT);
	$newList = $newList . "A$num_m.fa\t$list\n";
}
open(OUT, ">$outPath/contigfiles.list");
print OUT ("$newList");
close(OUT);

##### running prokka #####
mkdir("$outPath/Assembly_ann", 0755) || die "$!" if(!-e "$outPath/Assembly_ann");

my @list = `ls $outPath/Assembly`;
foreach my $list(@list) {
	chomp($list);
	my $fn = $list;
	$fn =~ s/\.\w+$//g;
	system("/opt/prokka/bin/prokka --outdir $outPath/Assembly_ann/$fn --prefix $fn $outPath/Assembly/$list --cpus $threads");
}
##### #####

##### copy ffn & gff files #####
mkdir("$outPath/FFN", 0755) || die "$!" if(!-e "$outPath/FFN");
mkdir("$outPath/GFF", 0755) || die "$!" if(!-e "$outPath/GFF");

system("cp $outPath/Assembly_ann/*/*.ffn $outPath/FFN/");
system("cp $outPath/Assembly_ann/*/*.gff $outPath/GFF/");
##### #####

