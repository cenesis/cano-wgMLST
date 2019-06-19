## cano-wgMLST files ##

1. Open Virtualization Format file: [cano-wgMLST.ova.7z](http://baccompare.imst.nsysu.edu.tw/download/cano-wgMLST.ova.7z) (5.2 GB)
2. The source file: [cano-wgMLST.7z](http://baccompare.imst.nsysu.edu.tw/download/cano-wgMLST.7z) (3 MB)
3. README file: [README.md](http://baccompare.imst.nsysu.edu.tw/download/README.md) (5 KB)
4. Results of _C. jejuni_: [C.jejuni.7z](http://baccompare.imst.nsysu.edu.tw/download/C.jejuni.7z) (59 MB)

## How to run ##

1. run with VirtualBox or VMware

	Use VirtualBox (<https://www.virtualbox.org>) or VMware (<http://www.vmware.com>) to import "cano-wgMLST.ova".
	The root password is 'rootwgmlst', the normal user (wgmlst) password is 'wgmlst'.

		$ 7za x cano-wgMLST.ova.7z

	or

2. install cano-wgMLST (cano-wgMLST.7z)

	To install all related tools, and use Perl and Python scripts to run.
		
	###Installation###

		Use CentOS 7.4.1708 as OS to show how to install.

		$ 7za x cano-wgMLST.7z
		$ cd cano-wgMLST

		* set PATH
			setenv PATH /YOUR_DIR/cano-wgMLST:/YOUR_DIR/cano-wgMLST/bin:$PATH

	###Software required###

		Following is a list of the softwares that you need:

		* Perl v5.16.3 (https://www.perl.org)
		* Python 2.7.5 (https://www.python.org)
		* Prokka v1.13 (https://github.com/tseemann/prokka)
		* Roary v3.12.0 (https://github.com/sanger-pathogens/Roary)
		* R v3.5.0 (https://www.r-project.org)
		* Fastx-Toolkit v0.0.14 (http://hannonlab.cshl.edu/fastx_toolkit/download.html)
		* PHYLIP v3.696 (http://evolution.genetics.washington.edu/phylip.html)
		* pip (https://pypi.org/project/pip/)
			$ wget https://bootstrap.pypa.io/get-pip.py
			$ python get-pip.py
		* ETE3-Toolkit (http://etetoolkit.org)
			$ sudo pip install ete3
		* scikit-learn (http://scikit-learn.org/stable/)
			$ sudo pip install -U scikit-learn
		* pdftk
			$ sudo yum localinstall https://www.linuxglobal.com/static/blog/pdftk-2.02-1.el7.x86_64.rpm
		* xvfb-run
			$ sudo yum install xorg-x11-server-Xvfb
		* IPython
			$ sudo yum install python-devel.x86_64 python27-python-devel.x86_64 gcc
			$ sudo pip install ipython


## Command details ##
### Run sample_run.sh ###
	
```
$ sample_run.sh
	Usage:
		sample_run.sh input_contigs_dir output_dir number_of_threads
	Arguments:
		Input contig directory [String]
		Output directory [String]
		Number of threads [Integer]
	Example: 
		Perform sample run by using 12 threads
		sample_run.sh ./sample_run/contigFiles ./sample_run 12
```


### Run 1.Contig-Annotator.pl ###
```
$ 1.Contig-Annotator.pl
	Usage:
		1.Contig-Annotator.pl -i input_dir -o output_dir [-p number_of_threads]
	Arguments:
		-i  Input contig directory [String]
		-o  Output annotated contig directory [String]
		-p  Number of threads [Integer]  Optional
		  default = 1
	Example: 
		Generate annotated contig files using 12 threads
		1.Contig-Annotator.pl -i ./contigFiles -o ./1.contigAnn -p 12
```


### Run 2.PGAdb-Builder.pl ###
```
$ 2.PGAdb-Builder.pl 
	Usage:
		2.PGAdb-Builder.pl -i input_dir -o output_dir [-m min_identity] [-p number_of_threads]
	Arguments:
		-i  Input cintig annotated directory [String]
		-o  Output PGAdb directory [String]
		-m  Minimum percentage identity for blastp [Integer]  Optional
		  default = 95
		-p  Number of threads [Integer]  Optional
		  default = 1
	Example: 
		Generate PGAdb using 12 threads
		2.PGAdb-Builder.pl -i ./1.contigAnn -o ./2.PGAdb -p 12
```


### Run 3.wgMLST-Profiler.pl ###
```
$ 3.wgMLST-Profiler.pl 
	Usage:
	3.wgMLST-Profiler.pl -i input_dir -d profileDB_dir -o output_dir [-m align_identity] [-n align_coverage] [-p number_of_threads] [-s scheme_file]
	Arguments:
		-i  Input contig directory [String]
		-d  Input PGAdb directory [String]
		-o  Output wgProfiles directory [String]
		-m  Percentage of aligned identity for blastn [Integer]  Optional
		  default = 90
		-n  Aligned coverage for blastn [Real]  Optional
		  default = 0.9
		-p  Number of threads [Integer]  Optional
		  default = 1
		-s  Selected scheme file [String]  Optional
	Example:
		Generate whole genome profiles using 12 threads
		3.wgMLST-Profiler.pl -i ./contigFiles -d ./2.PGAdb -o ./3.wgProfiles -p 12
```


### Run 4.Dendro-Plotter.pl ###
```
$ 4.Dendro-Plotter.pl 
	Usage:
		4.Dendro-Plotter.pl -i input_dir -o output_dir [-s scheme_file]
	Arguments:
		-i  Input wgProfiles directory [String]
		-o  Output dendrogram directory [String]
		-s  Selected scheme file [String]  Optional
		  default = 'input_dir/scheme/core.scheme'
	Example:
		Dendro-Plotter with core scheme
		4.Dendro-Plotter.pl -i ./3.wgProfiles -o ./4.DendroPlot
```


### Run 5.Loci-Extractor.pl ###
```
$ 5.Loci-Extractor.pl
	Usage:
		5.Loci-Extractor.pl -i input_dir -o output_dir -s scheme_file [-t top_threshold]
	Arguments:
		-i  Input dendrogram directory [String]
		-o  Output extracted loci directory [String]
		-s  Selected scheme file [String]  Optional
		  default = 'input_dir/core.scheme'
		-t  Top threshold [Integer]  Optional
		  default = 5
	Example: 
		Loci-Extractor
		5.Loci-Extractor.pl -i ./4.DendroPlot -o ./5.ExtractedLoci 
```

```
Date: 2019/06/17
Author:	Yen-Yi Liu (current788@gmail.com),
	Ji-Wei Lin (jwlin@imst.nsysu.edu.tw),
	Chih-Chieh Chen (chieh@imst.nsysu.edu.tw)
```
