#!/usr/bin/perl

=head1 NAME

 get_happy_matrix.pl
 This script built the matrix of ocurrence of each marker in each aliquot from a Happy Mapping or Radiation Hybrid experiment.
 As input, a preprocessed fasta file with aliquot tags at the beginning of the name is required.

=cut

=head1 AUTHORS

  Noe Fernandez-Pozo
  (nf232@cornell.edu)

=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

################################################################################
#                               LOAD & CHECKING                                #
################################################################################


if (!$ARGV[0]) {
	help_msg();
	exit;
}

# V0.0.2 includes rh_tsp_map input matrix as a result and a statistics file
# V0.0.3 is ready to accept any number of aliquots and is independent of the read name
# v0.0.4 memory usage was improved fixing duplicated variables that were referenced
my $VERSION = "0.0.4";

# ------------------------------------------------------------------------------ store options in hash
my %Options=get_options_and_params();

# ------------------------------------------------------------------------------ call to the help
if (($Options{help})||(!%Options)) {
	help_msg();
	exit;
}
# ------------------------------------------------------------------------------ ask the version
if ($Options{version}) {
	version_msg();
	exit;
}
# ------------------------------------------------------------------------------ check everything before start
check_options(\%Options);
# ------------------------------------------------------------------------------ saving and checking input subroutines

# save options and parameters in a hash
sub get_options_and_params {
	my %h=();
	# set default values
	$h{length} = 69;
	$h{reads} = 3;
	$h{aliquots} = 3;
	$h{percentage} = 20;
	
	GetOptions( \%h,
		"aliquots|a=i",
		"fasta|f=s",
		"tags|t=s",
		"help|h",
		"length|l=i",
		"percentage|p=i",
		"reads|r=i",
		"verbose",
		"version|v"
	);
	return %h;
}

# check if the options and parameters are right and sufficient for the execution
sub check_options {
	my $options = shift;
	
	# -------------------------------------------- check mandatory options
	if (!$$options{fasta}) {
		print STDERR "\n---> ERROR: a fasta file must be provided <---\n\n";
		
		help_msg();
		exit;
	}
	
	if (!$$options{tags}) {
		print STDERR "\n---> ERROR: a list of tags must be provided <---\n\n";
		
		help_msg();
		exit;
	}
	
	# -------------------------------------------- check files
	check_file_exist($$options{fasta},"---> File ".$$options{fasta}." was not found. A file with the marker reads is required <---");
	check_file_exist($$options{tags},"---> File ".$$options{tags}." was not found. A file with tags used for each aliquot of the happy mapping is required:\n\nacag\nactg\nagac\n...\n");
	
	# -------------------------------------------- removing old files and making new directories
	if (-e "output/") {
		`rm -rf output`;
	}
	mkdir ("output");
	
	#-------------------------------------------- print option chosed
	if ($$options{verbose}) {
		version_msg();
		
		print STDOUT "\n-------------- PARAMETERS AND OPTION CHOSEN --------------\n";
		foreach my $opt (keys %$options) {
			my $sp = (10 - length($opt));
			print STDOUT "${opt}".(" "x$sp)."  ==>  $$options{$opt}\n";
		}
		print STDOUT "----------------------------------------------------------\n";
	}
	
	#-------------------------------------------- print option chosed in a stats file
	open (my $stats_fh, ">output/get_happy_matrix_stats.txt") || die ("\nERROR: the file 'output/get_happy_matrix_stats.txt' could not be created\n");
	
	print $stats_fh "\n$0 v$VERSION\n";
	
	print $stats_fh "\n-------------- PARAMETERS AND OPTION CHOSEN --------------\n";
	foreach my $opt (keys %$options) {
		my $sp = (10 - length($opt));
		print $stats_fh "${opt}".(" "x$sp)."  ==>  $$options{$opt}\n";
	}
	print $stats_fh "----------------------------------------------------------\n";
	close($stats_fh);
	
	if ($$options{verbose}) {
		print STDOUT "--- YOUR PARAMETERS AND OPTIONS WERE CORRECTLY CHECKED ---\n";
	}
}

sub check_file_exist {
	my $file = shift;
	my $die_msg = shift;
	
	if (!-e $file) {
		print STDERR "\n$die_msg\n";
		
		help_msg();
		exit;
	}
}

sub version_msg {
	print STDOUT "\n$0 v$VERSION\n"
}

sub help_msg {

    print <<END;

Usage: $0 -f <file.fasta> -t <tags_file.txt> [other options or parameters]

  *** MANDATORY OPTIONS ****************************************************************
  
  -fasta|f <file>           A fasta file with the processed happy mapping reads labeled 
                            with their aliquot tag at the beginning of the name.
  -tags|t <file>            A text file with the happy mapping tags.
 
  *** MORE OPTIONS AND PARAMETERS ******************************************************

  -aliquots|a <integer>     Minumum number of aliquots by each marker to be included in
                            the matrix (default=3).
  -length|l <integer>       Minimum length of the markers (default=69bp).
  -percentage|p <integer>   Minimum percentage of reads by aliquot in relation to the average
                            to include the aliquot in the matrix (default=20).
  -reads|r <integer>        Minumum number of reads by each marker to be included in
                            the matrix (default=3).
  
  -verbose                  Verbose mode
  -version|v                Shows program version
  -help|h                   Shows this help

END

}


################################################################################
#                                 SUBROUTINES                                  #
################################################################################

# loading tags from Happy Tags file in a hash, saving their position in the matrix in the value of the hash
sub load_tags{
	
	my $tags_file = shift;
	my %tags_h;
	my $counter = 0;
	
	open (my $tag_fh, $tags_file) || die ("\nERROR: the file '$tags_file' could not be found\n");
	while (my $line = <$tag_fh>) {
		chomp($line);
		$tags_h{uc($line)} = $counter;
		$counter++;
	}
	return (\%tags_h);
}

# read a fasta file and save all the info about the markers and aliquots
sub read_fasta {
	
	my $opts = shift;
	my $tags = shift;
	
	my $counter = 0;
	my $tag_seq = '';
	my $name = '';
	my %markers; # save markers and their ocurrence in each aliquot
	my %marker_name; # save markers names and their sequences
	my @tag_num; # array with one position by aliquot
	my $tag_position = 0; # variable to save the position of each aliquot
	my $seq = '';
	my %reads_by_aliquot;
	my $reads_with_n_counter = 0;
	my $max_reads_with_n = 0;
	my %discarded_aliquots;
	
	if ($$opts{verbose}) {
		print "Reads  Markers\n";
	}
	
	my $input_fasta = Bio::SeqIO->new(-file => $$opts{fasta}, -format => 'fasta');
	
	while(my $seq_obj = $input_fasta->next_seq()) {
		
		# every read is counted
		$counter++;
		my $seq_name = $seq_obj->id();
		my $seq_fasta = $seq_obj->seq();
	
		if ($seq_name =~ m/^([a-t]{4}\d*)_(.+)/i) {
			$tag_seq = $1;
			$name = $2;
			$reads_by_aliquot{$tag_seq}++;
		} else {
			print "wrong sequence name: $seq_name\n";
		}
	
		if ($seq_fasta =~ m/N+/i) {
			$reads_with_n_counter++;
			my $n_num = $seq_fasta =~ tr/N/N/;

			if ($n_num > $max_reads_with_n) {
				$max_reads_with_n = $n_num;
				if ($$opts{verbose}) {
					print STDOUT "\nSample with $n_num Ns: $seq_fasta\n\n";
				}
			}
			next; # reads containing Ns are discarded
		}
		if (length($seq_fasta) >= $$opts{length}) {
			$seq = substr($seq_fasta,0,($$opts{length}-1));
		} else {
			print STDOUT "$seq_name was shorter than ".$$opts{length}." bp\n";
			next;
		}
		
		# ask for the aliquot position in the matrix (column number)
		$tag_position = int($$tags{uc($tag_seq)});
		
		# save markers in hash and count their ocurrence in each aliquot
		if (!$markers{$seq}) {
			my @tag_num=(0) x scalar(keys(%$tags));
			$markers{$seq} = \@tag_num;
			$markers{$seq}[$tag_position] = 1;
			$marker_name{$seq}=$seq_name;
		} else {
			$markers{$seq}[$tag_position]++;
		}
		# print "$tag: pos: $tag_position --> @{$markers{$seq}}\n";
		
		if ($$opts{verbose}) {
			if ($counter % 100000 == 0) {
				print "$counter ".scalar(keys %markers)."\n";
			}
		}
	}
	my $aliquot_read_mean = sprintf("%.2f",($counter/scalar(keys %reads_by_aliquot)));
	
	if ($$opts{verbose}) {
		print STDOUT "$counter ".scalar(keys %markers)."\n";
		print STDOUT "\n$reads_with_n_counter reads containing Ns were discarded\n";
		print STDOUT "Maximum number of Ns in a read: $max_reads_with_n\n";
		print STDOUT "\nAverage reads found per aliquot before filtering: $aliquot_read_mean\n";
	}
	#-------------------------------------------- opening stats file to append
	open (my $stats_fh, ">>output/get_happy_matrix_stats.txt") || die ("\nERROR: the file 'output/get_happy_matrix_stats.txt' could not be opened\n");
	print $stats_fh "\n$reads_with_n_counter reads containing Ns were discarded\n";
	print $stats_fh "Maximum number of Ns in a read: $max_reads_with_n\n";
	print $stats_fh "\nAverage reads found per aliquot before filtering: $aliquot_read_mean\n";
	#--------------------------------------------
	
	print STDOUT "\nReads distribution by aliquot before filtering\n";
	print STDOUT "\nAliquot\tReads\n";
	print $stats_fh "\nReads distribution by aliquot before filtering\n";
	print $stats_fh "\nAliquot\tReads\n";
	
	
	foreach my $aliq (keys %reads_by_aliquot) {
		
		my $aliq_percent = sprintf("%.2f",($reads_by_aliquot{$aliq}*100/$aliquot_read_mean));
		
		if ($aliq_percent < $$opts{percentage}) {
			if ($$opts{verbose}) {
				print STDOUT "$aliq\t$reads_by_aliquot{$aliq} --> $aliq_percent% of average, $aliq will be removed\n";
			}
			$discarded_aliquots{int($$tags{uc($aliq)})} = 1;
			
			#------------- appending to stats file discarded aliquots
			print $stats_fh "$aliq\t$reads_by_aliquot{$aliq} --> $aliq_percent% of average, $aliq will be removed\n";
		} else {
			print STDOUT "$aliq\t$reads_by_aliquot{$aliq}\n";
			print $stats_fh "$aliq\t$reads_by_aliquot{$aliq}\n";
		}
	}
	
	close($stats_fh);
	return (\%markers,\%marker_name,\%discarded_aliquots);
}


sub print_matrix {
	
	my $opts = shift;
	my $tags = shift;
	my $markers = shift;
	my $marker_name = shift;
	my $discarded_aliquots = shift;
	
	my $matrix_line = '';
	my @aliquot_reads = (0)x49;
	my $useful_markers = 0;
	my $aliquots_num = 0;
	my $matrix_head='locus_name';
	my @matrix_text;
	my @rh_matrix_text;
	
	my $m_num = sprintf("%07d", 0);
	open (my $markers_fh, ">output/markers_seqs.txt") || die ("\nERROR: the file 'output/markers_seqs.txt' could not be created\n");
	
	# each marker
	foreach my $marker (keys %$markers) {
		
		# each aliquot
		foreach my $i (0 .. $#{$$markers{$marker}}) {
			
			if (!$$discarded_aliquots{$i}) {
				if ($$markers{$marker}[$i] >= $$opts{reads}) {
					$$markers{$marker}[$i] = 1;
					$aliquots_num++;
				} else {
					$$markers{$marker}[$i] = 0;
				}
				$matrix_line = "${matrix_line}$$markers{$marker}[$i]";
			}
		}
		if ($aliquots_num >= $$opts{aliquots}) {
			$m_num++;
			print $markers_fh "m$m_num\t$$marker_name{$marker}\t$marker\n";
			$useful_markers++;
			
			# save data in Carthagene format
			push(@matrix_text,"*m$m_num\t$matrix_line");
			# save data in rh_tsp_map format
			my $id = $m_num;
			$id =~ s/^0+//;
			push(@rh_matrix_text,"$id\tm$m_num\t$matrix_line");
		}
		$matrix_line ='';
		$aliquots_num =0;
	}
	# ------------------------------------------------------------------------------------------------ Lets print the matrix
	open (my $matrix_fh, ">output/carthagene_input.txt") || die ("\nERROR: the file 'output/carthagene_input.txt' could not be created\n");
	open (my $rh_matrix_fh, ">output/rh_tsp_map_input.txt") || die ("\nERROR: the file 'output/rh_tsp_map_input.txt' could not be created\n");
	
	# print Carthagene input
	print $matrix_fh "data type radiated hybrid\n";
	print $matrix_fh "".(scalar(keys(%$tags)) - scalar(keys(%$discarded_aliquots)))." $useful_markers 0 0 0=A 1=H\n";
	print $matrix_fh join("\n", @matrix_text);
	
	# print rh_tsp_map input
	print $rh_matrix_fh "Id\tName\tRHvector\n";
	print $rh_matrix_fh join("\n", @rh_matrix_text);
	
	# ------------------------------------------------------------------------------------------------ opening stats file to append
	open (my $stats_fh, ">>output/get_happy_matrix_stats.txt") || die ("\nERROR: the file 'output/get_happy_matrix_stats.txt' could not be opened\n");
	
	print STDOUT "\n$useful_markers useful markers were selected after applying all filters:\n";
	print STDOUT "-->  Selected markers appeared at least in ".$$opts{aliquots}." aliquots\n";
	print STDOUT "-->  Selected markers were represented at least by ".$$opts{reads}." reads in these aliquots\n";
	print STDOUT "-->  Aliquots under ".$$opts{percentage}."% of the average of reads/aliquot were removed\n";
	print STDOUT "-->  Reads shorter than ".$$opts{length}." bp were discarded\n";
	print STDOUT "-->  Reads containing Ns were discarded\n";

	print $stats_fh "\n$useful_markers useful markers were selected after applying all filters:\n";
	print $stats_fh "-->  Selected markers appeared at least in ".$$opts{aliquots}." aliquots\n";
	print $stats_fh "-->  Selected markers were represented at least by ".$$opts{reads}." reads in these aliquots\n";
	print $stats_fh "-->  Aliquots under ".$$opts{percentage}."% of the average of reads/aliquot were removed\n";
	print $stats_fh "-->  Reads shorter than ".$$opts{length}." bp were discarded\n";
	print $stats_fh "-->  Reads containing Ns were discarded\n";

	close($markers_fh);
	close($matrix_fh);
	close($stats_fh);
}

################################################################################
#                                    MAIN                                      #
################################################################################

# ------------------------------------------------------------------------------ loading tags
my ($tags, $file_head) = load_tags($Options{tags});

# ------------------------------------------------------------------------------ read fasta
my ($markers, $markers_names, $bad_aliquots) = read_fasta(\%Options, $tags);

# ------------------------------------------------------------------------------ get the matrix
print_matrix(\%Options, $tags, $markers, $markers_names, $bad_aliquots);


