use strict;
use Cwd;
use FindBin;
my $curr_dir = getcwd;
my $myBin = $FindBin::Bin;
my $THREAD = 4;

my $SETTING_FILE = "setting";
my $HMMSCAN = "hmmscan";
my $PRODIGAL = "prodigal";
my $HMMSCANANN = "python $myBin/Hmmscan2Ann.py";
my $COGDB = "$myBin/data/COG_all.hmm";

checkProgram();

my $E_VALUE_CUTOFF = 1e-10;
my $COGLIST = "$myBin/COGlist.txt";
my $COGCATEGORY = "$myBin/COGcategory.txt";
my $COGLIST2 = "$myBin/COG.list";

if (!(defined($ARGV[0]) && defined($ARGV[1])))
{
	print "Usage: perl map_COG.pl (list file) [(output file)]\n";
	exit; 
}

my $cog_f = $ARGV[0];
my $out_f = $ARGV[1];

mapCOG($cog_f, $out_f);

sub mapCOG
{
	my $cog_f = $_[0];
	my $out_f = $_[1];

	my %category_hash;
	my %COGCategory_hash;

	my @samplearr;
	my $sample;

	my %COGDBhash;

	my $line;
	my $cmd;
	my $tmpline;
	my %tmphash;
	my @arr;
	my @tmparr;
	my $arrnum;
	my $i;
	my $j;
	my $k;
	my $l;
	my $p;

	# Read COGLIST file
	open(FILE, "<$COGLIST");
	while(defined($line = <FILE>))
	{
		chomp($line);
		if ($line ne "")
		{
			@arr = split(/\t/, $line);
			$COGDBhash{$arr[0]} = $arr[1];
		}
	}
	close(FILE);

	# Read COG category file
	open(FILE, "<$COGCATEGORY");
	while(defined($line = <FILE>))
	{
		chomp($line);
		@arr = split(/\t/, $line);
		$category_hash{$arr[0]} = $arr[1];
	}
	close(FILE);

	# Read COG list2 file (category file)
	open(FILE, "<$COGLIST2");
	while(defined($line = <FILE>))
	{
		chomp($line);
		if ($line =~ /\[([A-Z]+)\]\t([A-Za-z0-9]+)/)
		{
			$COGCategory_hash{$2} = $1;
		}
	}
	close(FILE);

	# Read list
	open(FILE, "<$cog_f");
	while(defined($line = <FILE>))
	{
		chomp($line);
		if ($line ne "")
		{
			push(@samplearr, $line);
		}
	}
	close(FILE);

	# map the aa sequences against COG database
	foreach $sample (@samplearr)
	{
		if (-e $sample)
		{
			print "Processing file $sample\n";
			if (!(-e "$sample.COG_annotation"))
			{
				$i = checkFormat($sample);
				if ($i == 1)
				{
					print "--Predicting genes from $sample\n";
					$cmd = "$PRODIGAL -a $sample.faa -i $sample -p meta 1>/dev/null 2>/dev/null";
					system($cmd);
					$j = "$sample.faa";
				}
				else
				{
					$j = $sample;
				}
				print "--Annotating genes in terms of COGs\n";
				$cmd = "$HMMSCAN --tblout $sample.hmmout -E $E_VALUE_CUTOFF --cpu $THREAD $COGDB $j 1>/dev/null 2>/dev/null";
				system($cmd);
				parseHMMout("$sample.hmmout", "$sample.COG_annotation");
				unlink("$sample.faa");
				unlink("$sample.hmmout");
			}
		}
		else
		{
			print "Cannot open file $sample.\n";
		}
	}

	my @outarr;
	my %filehash;
	my %tmphash;
	my $i = 0;
	my $k = 0;
	foreach $sample (@samplearr)
	{
		if (-e "$sample.COG_annotation")
		{
			$filehash{$sample} = $i;
			$i++;
			open(FILE, "<$sample.COG_annotation");
			while(defined($line = <FILE>))
			{
				@arr = split(/\t/, $line);
				if ($arr[1] =~ /.(COG[0-9]+)./)
				{
					$j = $COGCategory_hash{$1};
					@tmparr = split(//, $j);
					foreach $j (@tmparr)
					{
						if (!(exists $tmphash{$j}))
						{
							$tmphash{$j} = $k;
							$k++;
						}
						if (defined($outarr[$filehash{$sample}][$tmphash{$1}]))
						{
							$outarr[$filehash{$sample}][$tmphash{$j}]++;
						}
						else
						{
							$outarr[$filehash{$sample}][$tmphash{$j}] = 1;
						}
					}
				}
			}
			close(FILE);
		}
	}

	open(OUT, ">$out_f");
	foreach $sample (@samplearr)
	{
		print OUT "\t$sample";
	}
	print OUT "\n";

	foreach $i (sort keys %tmphash)
	{
		print OUT "$i\t$category_hash{$i}";
		foreach $sample (@samplearr)
		{
			if (defined($outarr[$filehash{$sample}][$tmphash{$i}]))
			{
				print OUT "\t$outarr[$filehash{$sample}][$tmphash{$i}]";
			}
			else
			{
				print OUT "\t0";
			}
		}
		print OUT "\n";
	}
}

sub parseHMMout
{
	my $hmmout_f = $_[0];
	my $out_f = $_[1];

	my $line;
	my @arr;
	my %temphash;
	my $curr = "";

	open(HMMFILE, "<$hmmout_f") || die "Cannot open HMM output file $hmmout_f. Please check if your HMMER is installed correctly.\n";
	open(HMMOUT, ">$out_f");
	while(defined($line = <HMMFILE>))
	{
		if ($line =~ /^\#/)
		{
			next;
		}
		@arr = split(/[ ]+/, $line);
		if ($arr[2] ne $curr)
		{
			$curr = $arr[2];
			print HMMOUT "$arr[2]\t$arr[0]\n";
		}
	}
	close(HMMFILE);
	close(HMMOUT);
}

# (Read and Check program directory)
sub checkProgram
{
	my $localline;
	my $tmpstr;
	my $tmpname = "tmp_" . time();

	if (!(-e $COGDB))
	{
		print "Please download the COG hmm files from eggNOG website and concatenate them into COG_all.hmm as instructed in the tutorial PDF file\nPlease place them under \'data\' folder at the COGmapper folder.\n";
		exit;
	}
	if (!(-e "$COGDB.h3i"))
	{
		print "Please compress the COG_all.hmm file";
		exit;
	}

	open(FILE, "<$myBin\/$SETTING_FILE");
	while(defined($localline = <FILE>))
	{
		chomp($localline);
		if ($localline =~ /\[([A-Za-z0-9_]+)\] ([A-Za-z0-9._\(\)\[\]\{\}\|\$\!\=\-\+\\\/]+)/)
		{
			if ($1 eq "Prodigal")
			{
				if (-d $2 && -e "$2\/$PRODIGAL")
				{
					$PRODIGAL = $2 . "\/" . $PRODIGAL;
				}
			}
			elsif ($1 eq "HMMER3")
			{
				if (-d $2 && -e "$2\/$HMMSCAN")
				{
					$HMMSCAN = $2 . "\/" . $HMMSCAN;
				}
			}
		}
	}
	close(FILE);

	# Check program
	# Prodigal
	$localline = "$PRODIGAL 1>/dev/null 2>$tmpname";
	system($localline);
	$tmpstr = "";
	open(FILE, "<$tmpname");
	while(<FILE>)
	{
		$tmpstr .= $_;
	}
	if ($tmpstr !~ /PRODIGAL/)
	{
		print "Cannot run Prodigal. Please indicate the file directory in \'$SETTING_FILE\' file.\n";
		exit;
	}

	# HMMER3
	$localline = "$HMMSCAN 1>$tmpname 2>/dev/null";
	system($localline);
	$tmpstr = "";
	open(FILE, "<$tmpname");
	while(<FILE>)
	{
		$tmpstr .= $_;
	}
	if ($tmpstr !~ /hmmscan/)
	{
		print "Cannot run HMMER3. Please indicate the file directory in \'$SETTING_FILE\' file.\n";
		exit;
	}

	unlink("$tmpname");
}

sub checkFormat
{
	my $f = $_[0];
	my $myline; 
	my $start = 0;
	my $m;
	my $total = 0;
	my $ATCG = 0;

	open(MYFILE, "<$f");
	while(defined($myline = <MYFILE>))
	{
		chomp($myline);
		if ($myline eq "")
		{
			next;
		}
		if ($start == 0)
		{
			if ($myline =~ /^>/)
			{
				$start = 1;
			}
			else
			{
				print "$f does not seem to be in FASTA format. Program stop.\n";
				exit;
			}
		}
		else
		{
			if ($myline !~ /^>/)
			{
				$m = length($myline);
				$total += $m;
				$myline =~ s/[ATCGatcg]//g;
				$m = $m - length($myline);
				$ATCG += $m;
			}
		}
	}
	close(MYFILE);

	$m = $ATCG / $total;
	if ($m > 0.8)
	{
		return 1;
	}
	else
	{
		return 2;
	}
}

