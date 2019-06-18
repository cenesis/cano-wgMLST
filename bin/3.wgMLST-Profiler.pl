#!/usr/bin/perl

use File::Basename;
my $programName = basename($0);

use Getopt::Long;
my $usage = "
Usage:   $programName -i input_dir -d profileDB_dir -o output_dir [-m align_identity] [-n align_coverage] [-p number_of_threads] [-s scheme_file]

Arguments: -i  Input contig directory [String]
           -d  Input PGAdb directory [String]
           -o  Output wgProfiles directory [String]
           -m  Percentage of aligned identity for blastn [Integer]  Optional
             default = 90
           -n  Aligned coverage for blastn [Real]  Optional
             default = 0.9
           -p  Number of threads [Integer]  Optional
             default = 1
           -s  Selected scheme file [String]  Optional

Example: Generate whole genome profiles using 12 threads
         $programName -i ./contigFiles -d ./2.PGAdb -o ./3.wgProfiles -p 12\n\n\n";


my $inPath = undef;
my $DB_path = undef;
my $outPath = undef;
my $pident_cut = undef;
my $aligcov_cut = undef;
my $threads = undef;
my $userScheme = undef;

die $usage unless GetOptions(
		'i|input_dir=s'          => \$inPath,
		'd|db_dir=s'             => \$DB_path,
        'o|output_dir=s'         => \$outPath,
        'm|pident=i'             => \$pident_cut,
        'n|aligcov=f'            => \$aligcov_cut,
        'p|threads=i'            => \$threads,
        's|scheme_file=s'        => \$userScheme)
	&& defined $inPath 
	&& defined $DB_path
	&& defined $outPath
#	&& defined $pident_cut
#	&& defined $aligcov_cut
#	&& defined $threads
#	&& defined $userScheme
	&& @ARGV == 0;

$pident_cut = 90 if($pident_cut == undef);
$aligcov_cut = 0.9 if($aligcov_cut == undef);
$threads = 1 if($threads == undef);

opendir(DIR, "$inPath") || die "\nError in opening contigs directory $inPath\n\n";
my @contig_files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;
closedir(DIR);
my $contig_num = @contig_files;
die "\nPlease input at least five genome contig files.\n\n" if($contig_num < 5);

open(DB,"<","$DB_path/locusfiles.list") || die "\nError in opening profileDB directory $DB_path\n\n";
close(DB);

mkdir("$outPath", 0755) || die "$!" if(!-e "$outPath");

#exit();


##### reformat contig files #####

mkdir("$outPath/Assembly", 0755) || die "$!" if(!-e "$outPath/Assembly");
my @list = `ls $inPath`;
my $m = 0;
my $newList = "";
my $matrixLabel = "#";
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

	my $isolateName = $list;
	$isolateName =~ s/\.\w+$//g;
	$matrixLabel = $matrixLabel . "\t$isolateName";
}
$matrixLabel = $matrixLabel . "\n";

open(OUT, ">$outPath/contigfiles.list");
print OUT ("$newList");
close(OUT);

##########

##### copy locusmapping #####

system("cp -rf $DB_path/locusmapping $outPath/") if(-e "$DB_path/locusmapping");

##########

##### copy scheme #####

system("cp -rf $userScheme $outPath/") if(length($userScheme) != 0);


##########

##### read userScheme file #####
my %userScheme = ();

if(length($userScheme) != 0) {
	my @scheme = `cat $userScheme`;
	foreach my $scheme(@scheme) {
		chomp($scheme);
		$userScheme{$scheme} = 1;
	}
}
##########

##### profiling #####

system("rm -rf $outPath/panDB_all.fa") if(-e "$outPath/panDB_all.fa");

if(length($userScheme) != 0) {
	my $panDB = "";

	my $flag = 0;
	system("cat $DB_path/locusfiles/* > $outPath/panDB_all_tmp.fa");
	open(PANDB, "$outPath/panDB_all_tmp.fa") or die("Could not open panDB_all_tmp.fa file.");
	foreach my $line(<PANDB>){
		if($line =~ /^>/) {
			my @temp = split(/\:\:/, $line);
			$temp[0] =~ s/^>//g;
			if($userScheme{$temp[0]} == 1) {
				$panDB = $panDB . $line;
				$flag = 1;
			} else {
				$flag = 0;
			}
		} else {
			if($flag == 1) {
				$panDB = $panDB . $line;
				$flag = 0;
			}
		}
	}
	close(PANDB);
	system("rm -rf $outPath/panDB_all_tmp.fa") if(-e "$outPath/panDB_all_tmp.fa");

	open(OUT, ">$outPath/panDB_all.fa");
	print OUT ("$panDB");
	close(OUT);
	$panDB = "";

} else {
	system("cat $DB_path/locusfiles/* > $outPath/panDB_all.fa");
}

  ##### read panDB_all.fa sequence length #####

my %panLen = ();
my $locusName = "";
my $alleleNum = "";

open(SEQ, "$outPath/panDB_all.fa") or die("Could not open panDB_all.fa file.");
foreach my $line(<SEQ>) {
	chomp($line);
	if($line =~ /^>/) {
		my @temp = split(/\:\:/, $line);
		$temp[0] =~ s/^>//g;
		$locusName = $temp[0];
		$alleleNum = $temp[1];
	} else {
		$panLen{$locusName}{$alleleNum} = length($line);
	}
}
close(SEQ);

  ##########

my @locusFiles = `ls $DB_path/locusfiles`;
my $locusNum = @locusFiles;
my $assemblyNum = $m;

##########


for(my $n=1;$n<=$assemblyNum;$n++) {

	my $asmFileName = (sprintf "A%08d.fa", $n);

	system("cp $outPath/Assembly/$asmFileName $outPath/AssemblyDB.fa");

	system("makeblastdb -in $outPath/AssemblyDB.fa -dbtype nucl -out $outPath/AssemblyDB");
	system("blastn -num_threads $threads -query $outPath/panDB_all.fa -db $outPath/AssemblyDB -out $outPath/profiling.bls -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'");


	##### read profiling.bls #####

	my %locusDB = ();
	my %alleleDB = ();
	my %qlenMax = ();

	open(BLS, "$outPath/profiling.bls") or die("Could not open profiling.bls file.");
	foreach my $line(<BLS>){
		chomp($line);
		my ($qseqid,$sseqid,$pident,$length,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore) = split(/\t/, $line);
		my @tmpQ = split(/\:\:/, $qseqid);
		my @tmpS = split(/\:\:/, $sseqid);
		my $qlen = $panLen{$tmpQ[0]}{$tmpQ[1]};
		my $aligcov = (($length-$gapopen)/$qlen);

		if($aligcov >= $aligcov_cut && $pident >= $pident_cut) {
			$locusDB{$tmpQ[0]} = 1;
		}

		if ($pident == 100 && $length == $qlen && $mismatch == 0 && $gapopen == 0 ) {
			if($length > $qlenMax{$tmpQ[0]}) {
				$alleleDB{$tmpQ[0]} = $tmpQ[1];
				$qlenMax{$tmpQ[0]} = $length;
			}
		}
	}
	close(BLS);

	##########




	##### print locus matrix #####

	my $locusMatrix = "";

	for(my $i=1;$i<=$locusNum;$i++) {
		my $locusID = (sprintf "SAL%07d", $i);
		next if(length($userScheme) != 0 && $userScheme{$locusID} != 1);
		$locusMatrix = $locusMatrix . $locusID;
		$locusDB{$locusID} = 0 if($locusDB{$locusID} eq "");
		$locusMatrix = $locusMatrix . "\t" . $locusDB{$locusID} . "\n";
	}
	open(OUT, ">$outPath/locusMatrix.$n");
	print OUT ("$locusMatrix");
	close(OUT);

	##########


	##### print alleles matrix #####

	my $alleleMatrix = "";

	for(my $i=1;$i<=$locusNum;$i++) {
		my $locusID = (sprintf "SAL%07d", $i);
		next if(length($userScheme) != 0 && $userScheme{$locusID} != 1);
		$alleleMatrix = $alleleMatrix . $locusID;
		$alleleDB{$locusID} = 0 if($alleleDB{$locusID} eq "");
		$alleleMatrix = $alleleMatrix . "\t" . $alleleDB{$locusID} . "\n";
	}

	open(OUT, ">$outPath/alleleMatrix.$n");
	print OUT ("$alleleMatrix");
	close(OUT);

	##########
}

##### marge all locusMatrix.xxx files to locusMatrix #####

my %locusHash = (); 

for(my $i=1;$i<=$assemblyNum;$i++) {
	my $fileName = "locusMatrix.$i";
	my @locus = `cat $outPath/$fileName`;
	foreach my $locus(@locus) {
		chomp($locus);
		my @temp = split(/\t/, $locus);
		$locusHash{$temp[0]} = $locusHash{$temp[0]} . $temp[1] . "\t";
	}
}

open(OUT, ">$outPath/locusMatrix");
for(my $i=1;$i<=$locusNum;$i++) {
	my $locusID = (sprintf "SAL%07d", $i);
	next if(length($userScheme) != 0 && $userScheme{$locusID} != 1);
	$locusHash{$locusID} =~ s/\t$//g;
	print OUT ("$locusID\t$locusHash{$locusID}\n");
}
close(OUT);

open(OUT, ">$outPath/locus_profile.xls");
print OUT ("$matrixLabel");
for(my $i=1;$i<=$locusNum;$i++) {
	my $locusID = (sprintf "SAL%07d", $i);
	next if(length($userScheme) != 0 && $userScheme{$locusID} != 1);
	$locusHash{$locusID} =~ s/\t$//g;
	print OUT ("$locusID\t$locusHash{$locusID}\n");
}
close(OUT);

##########


##### marge all alleleMatrix.xxx files to alleleMatrix #####

my %alleleHash = (); 

for(my $i=1;$i<=$assemblyNum;$i++) {
	my $fileName = "alleleMatrix.$i";
	my @alleles = `cat $outPath/$fileName`;
	foreach my $alleles(@alleles) {
		chomp($alleles);
		my @temp = split(/\t/, $alleles);
		$alleleHash{$temp[0]} = $alleleHash{$temp[0]} . $temp[1] . "\t";
	}
}

open(OUT, ">$outPath/alleleMatrix");
for(my $i=1;$i<=$locusNum;$i++) {
	my $locusID = (sprintf "SAL%07d", $i);
	next if(length($userScheme) != 0 && $userScheme{$locusID} != 1);
	$alleleHash{$locusID} =~ s/\t$//g;
	print OUT ("$locusID\t$alleleHash{$locusID}\n");
}
close(OUT);

open(OUT, ">$outPath/allele_profile.xls");
print OUT ("$matrixLabel");
for(my $i=1;$i<=$locusNum;$i++) {
	my $locusID = (sprintf "SAL%07d", $i);
	next if(length($userScheme) != 0 && $userScheme{$locusID} != 1);
	$alleleHash{$locusID} =~ s/\t$//g;
	print OUT ("$locusID\t$alleleHash{$locusID}\n");
}
close(OUT);

##########



##### copy scheme and run varianceThreshold #####

system("cp -r $DB_path/scheme $outPath/");

use FindBin;
my $programPath = $FindBin::Bin;

for(my $i=5;$i<=10;$i++) {
	my $var_th = $i/10;
	system("python $programPath/varianceThreshold.py $outPath/alleleMatrix $var_th > $outPath/scheme/variance_$i.scheme");
}



#################################################
