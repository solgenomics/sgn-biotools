#!/usr/bin/perl

=head1 NAME

 get_happy_matrix.pl
 This script built the matrix of ocurrence of each marker in each aliquot from a happy mapping experiment.
 As input, a fasta file with aliquot tags at the beginning of the name is required.

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

my $VERSION = "0.0.1";

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
	my $opt = shift;
	my %options = %$opt;
	
	if ($options{verbose}) {
		version_msg();
		
		print STDERR "\n-------------- PARAMETERS AND OPTION CHOSEN --------------\n";
		foreach my $opt (keys %options) {
			my $sp = (10 - length($opt));
			print STDERR "${opt}".(" "x$sp)."  ==>  $options{$opt}\n";
		}
		print STDERR "----------------------------------------------------------\n";
	}
	
	# -------------------------------------------- check mandatory options
	if ((!$options{fasta}) || (!$options{tags})) {
		print "\nERROR: a fasta file and a list of tags must be provided\n\n";
		
		help_msg();
		exit;
	}
	
	# -------------------------------------------- check files
	check_file_exist($options{fasta},"File ".$options{fasta}." was not found. A file with the marker reads is required");
	check_file_exist($options{tags},"File ".$options{tags}." was not found. A file with tags used for each aliquot of the happy mapping is required:\n\nacag\nactg\nagac\n...\n");
	
	if ($options{verbose}) {
		print STDERR "--- YOUR PARAMETERS AND OPTIONS WERE CORRECTLY CHECKED ---\n";
	}
	
	# -------------------------------------------- removing old files and making new directories
	if (-e "output/") {
		`rm -rf output`;
	}
	mkdir ("output");
}

sub check_file_exist {
	my $file = shift;
	my $die_msg = shift;
	
	if (!-e $file) {
		print "\n$die_msg\n";
		
		help_msg();
		exit;
	}
}

sub version_msg {
	print "$0 v$VERSION\n"
}

sub help_msg {

    print <<END;

Usage: $0 -f <file.fasta> -t <tags_file.txt> [other options or parameters]

  *** MANDATORY OPTIONS ****************************************************************
  
  -fastq|f <file>           A fasta file with the processed happy mapping reads labeled 
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
	# my $head = 'locus_name';
	open (my $tag_fh, $tags_file) || die ("\nERROR: the file '$tags_file' could not be found\n");
	while (my $line = <$tag_fh>) {
		chomp($line);
		# $head = "$head\t".uc($line);
		$tags_h{uc($line)} = $counter;
		$counter++;
	}
	return (\%tags_h);
	# return (\%tags_h,$head);
}

# read a fasta file and save all the info about the markers and aliquots
sub read_fasta {
	
	my $opt_h = shift;
	my %opts = %$opt_h;
	
	my $tag_h = shift;
	my %tags = %$tag_h;
	
	my $input_fasta = Bio::SeqIO->new(-file => $opts{fasta}, -format => 'fasta');
	
	my $counter = 0;
	my $tag = '';
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
	
	if ($opts{verbose}) {
		print "Reads  Markers\n";
	}
	
	while(my $seq_obj = $input_fasta->next_seq()) {
		
		# every read is counted
		$counter++;
		my $seq_name = $seq_obj->id();
		my $seq_fasta = $seq_obj->seq();
	
		if ($seq_name =~ m/^([a-t]{4})_D3B4KKQ1\:291\:D17NUACXX\:8\:(.+)/i) {
			$tag = $1;
			$name = $2;
			$reads_by_aliquot{$tag}++;
		} else {
			print "wrong sequence name: $seq_name\n";
		}
	
		if ($seq_fasta =~ m/N+/i) {
			$reads_with_n_counter++;
			my $n_num = $seq_fasta =~ tr/N/N/;

			if ($n_num > $max_reads_with_n) {
				$max_reads_with_n = $n_num;
				if ($opts{verbose}) {
					print "\nSample with $n_num Ns: $seq_fasta\n\n";
				}
			}
			next; # reads containing Ns are discarded
		}
		if (length($seq_fasta) >= $opts{length}) {
			$seq = substr($seq_fasta,0,($opts{length}-1));
		} else {
			print "D3B4KKQ1\:291\:D17NUACXX\:8\:$seq_name was shorter than ".$opts{length}." bp\n";
			next;
		}
		
		# ask for the aliquot position in the matrix (column number)
		$tag_position = int($tags{uc($tag)});
		
		# save markers in hash and count their ocurrence in each aliquot
		if (!$markers{$seq}) {
			my @tag_num=(0)x48;
			$markers{$seq} = \@tag_num;
			$markers{$seq}[$tag_position] = 1;
			$marker_name{$seq}=$seq_name;
		} else {
			$markers{$seq}[$tag_position]++;
		}
		# print "$tag: pos: $tag_position --> @{$markers{$seq}}\n";
		
		if ($opts{verbose}) {
			if ($counter % 100000 == 0) {
				print "$counter ".scalar(keys %markers)."\n";
			}
		}
	}
	my $aliquot_read_mean = sprintf("%.2f",($counter/scalar(keys %reads_by_aliquot)));
	
	if ($opts{verbose}) {
		print "$counter ".scalar(keys %markers)."\n";
		print "\n$reads_with_n_counter reads containing Ns were discarded\n";
		print "Maximum number of Ns in a read: $max_reads_with_n\n";
		print "\nAverage reads found per aliquot before any filter: $aliquot_read_mean\n";
	}
	
	foreach my $aliq (keys %reads_by_aliquot) {
		
		my $aliq_percent = sprintf("%.2f",($reads_by_aliquot{$aliq}*100/$aliquot_read_mean));
		
		
		if ($aliq_percent < $opts{percentage}) {
			if ($opts{verbose}) {
				print "$aliq: $reads_by_aliquot{$aliq} found, $aliq_percent% of average, $aliq will be removed\n";
			}
			$discarded_aliquots{int($tags{uc($aliq)})} = 1;
		}
	}
	
	return (\%markers,\%marker_name,\%discarded_aliquots);
}


sub print_matrix {
	
	my $opt_h = shift;
	my $tag_h = shift;
	my $markers_h = shift;
	my $markers_names_h = shift;
	my $bad_aliq_h = shift;
	
	my %opts = %$opt_h;
	my %tags = %$tag_h;
	my %markers = %$markers_h;
	my %marker_name = %$markers_names_h;
	my %discarded_aliquots = %$bad_aliq_h;
	
	my $m_num = sprintf("%07d", 0);
	open (my $markers_fh, ">output/markers_seqs.txt") || die ("\nERROR: the file 'output/markers_seqs.txt' could not be created\n");

	my $matrix_line = '';
	my @aliquot_reads = (0)x49;
	my $useful_markers = 0;
	my $aliquots_num = 0;
	my $matrix_head='locus_name';
	my @matrix_text;
	
	foreach my $aliq (sort{$tags{$a} <=> $tags{$b}} keys %tags) {
		my $tmp_key = $tags{$aliq};
		if (!$discarded_aliquots{$tmp_key}) {
			$matrix_head = "$matrix_head\t".uc($aliq);
		# } else {
		# 	print "$aliq at $tags{$aliq} was no printed in the head\n";
		}
	}
	push(@matrix_text,$matrix_head);
	# print $matrix_fh "$matrix_head\n";
	
	# each marker
	foreach my $marker (keys %markers) {
		# print "@{$markers{$marker}}\n";
		# $m_num++;
		# print $markers_fh "m$m_num\t$marker_name{$marker}\t$marker\n";
		
		# each aliquot
		foreach my $i (0 .. $#{$markers{$marker}}) {
			if (!$discarded_aliquots{$i}) {
				if ($markers{$marker}[$i] >= $opts{reads}) {
					$aliquot_reads[$i]++;
					$markers{$marker}[$i] = 1;
					$aliquots_num++;
				} else {
					$markers{$marker}[$i] = 0;
				}
				$matrix_line = "${matrix_line}$markers{$marker}[$i]";
			}
		}
		if ($aliquots_num >= $opts{aliquots}) {
			$m_num++;
			print $markers_fh "m$m_num\t$marker_name{$marker}\t$marker\n";
			$useful_markers++;
			# print $matrix_fh "m${m_num}${matrix_line}\n";
			push(@matrix_text,"*m${m_num}\t${matrix_line}");
		}
		$matrix_line ='';
		$aliquots_num =0;
	}
	# ------------------------------------------------------------------------------------------------ Lets print the matrix
	open (my $matrix_fh, ">output/hm_matrix.txt") || die ("\nERROR: the file 'output/hm_matrix.txt' could not be created\n");
	print $matrix_fh "data type radiated hybrid\n";
	print $matrix_fh "".(scalar(keys(%tags)) - scalar(keys(%discarded_aliquots)))." $useful_markers 0 0 0=A 1=H\n";
	# print $matrix_fh join("\n", @matrix_text);
	# ----------------------------------------------------------------------------------------------------------------------
	
	print "\n";
	my $index=0;
	foreach my $aliq (sort{$tags{$a} <=> $tags{$b}} keys %tags) {
		if (!$discarded_aliquots{$index}) {
			print "".($index+1)."\t$aliq\t$aliquot_reads[$index]\n";
			# print "".($index+1)."\t$tags_array[$index]\t$aliquot_reads[$index]\n";
		}
		$index++;
	}
	print "\nUseful markers: $useful_markers\n";

	close($markers_fh);
	close($matrix_fh);
	
}

################################################################################
#                                    MAIN                                      #
################################################################################

# ------------------------------------------------------------------------------ loading tags
my ($tags_hash, $file_head) = load_tags($Options{tags});
my %tags = %$tags_hash;

# ------------------------------------------------------------------------------ read fasta
my ($markers_h,$markers_names_h,$bad_aliquots_h) = read_fasta(\%Options, \%tags);
my %markers = %$markers_h;
my %markers_names = %$markers_names_h;
my %bad_aliquots = %$bad_aliquots_h;

# ------------------------------------------------------------------------------ get the matrix
print_matrix(\%Options,\%tags,\%markers,\%markers_names,\%bad_aliquots);


