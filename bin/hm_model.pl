#!/usr/bin/perl

=head1 NAME

 hm_model
 This script simulates a happy maping experiment.

=cut

=head1 AUTHORS

  Noe Fernandez-Pozo
  (nf232@cornell.edu)

=cut

use strict;
use warnings;

use Getopt::Long;
use Bio::SeqIO;
use Math::Random;

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
	$h{dna_copies} = 100;
	$h{extraction_size} = 40000;
	$h{genome_percent} = 0.7;
	$h{std_dev} = 10000;
	$h{reads_length} = 69;
	$h{total_rs} = 'C^TGCAG';
	
	GetOptions( \%h,
		"aliquot_tags|a=s",
		"dna_copies|d=i",
		"extraction_size|e=i",
		"fasta|f=s",
		"genome_percent|g=f",
		"std_dev|s=i",
		"reads_length|r=i",
		"total_rs|t=s",
		"help|h",
		"verbose",
		"version|v"
	);
	return %h;
}

# check if the options and parameters are right and sufficient for the execution
sub check_options {
	my $opt = shift;
	my %options = %$opt;
	
	# -------------------------------------------- check mandatory options
	if (!$options{fasta}) {
		print "\nERROR: a fasta file must be provided\n\n";
		help_msg();
		exit;
	}
	if (!$options{aliquot_tags}) {
		print "\nERROR: a file with Happy Mapping Tags must be provided\n\n";
		help_msg();
		exit;
	}
	
	# -------------------------------------------- check files
	check_file_exist($options{fasta},"File ".$options{fasta}." was not found. A file with the marker reads is required");
	check_file_exist($options{aliquot_tags},"File ".$options{aliquot_tags}." was not found. A file with the marker reads is required");
	
	
	if ($options{verbose}) {
		version_msg();
		
		print STDERR "\n-------------- PARAMETERS AND OPTION CHOSEN --------------\n";
		foreach my $opt (keys %options) {
			my $sp = (10 - length($opt));
			print STDERR "${opt}".(" "x$sp)."  ==>  $options{$opt}\n";
		}
		print STDERR "----------------------------------------------------------\n";
	}
	
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

Usage: $0 -f <file.fasta> [other options or parameters]

  *** MANDATORY OPTIONS ****************************************************************
  
  -fasta|f <file>               A fasta file with the genome under study
  -aliquot_tags|a <file>        A file with Happy Mapping tags, one by line

  *** MORE OPTIONS AND PARAMETERS ******************************************************

  -dna_copies|d <integer>       Number of copies of the genome (default=100)
  -extraction_size|e <integer>  Average size of DNA after extraction (default=40000bp)
  -std_dev|s <integer>          Standard deviation for DNA size (default=10000bp)

  -genome_percent|g <float>     Percentage of genome included by aliquot (default=0.7)
  -total_rs|t <string>          Restriction site for a total digestion (default='C^TGCAG')
  -reads_length|r <integer>     Final reads length (default=69)
                                
  -verbose                      Verbose mode
  -version|v                    Shows program version
  -help|h                       Shows this help

END

}

################################################################################
#                                 SUBROUTINES                                  #
################################################################################

# Lets simulate a DNA fragmentation during the extraction based in an average size and a standard deviation
sub dna_extraction {
	my $opts = shift;
	
	my @dna;
	my $dna_size=0;
	my $input_fastq = Bio::SeqIO->new(-file => $$opts{fasta}, -format => 'fasta');
	
	# foreach read in the genome fasta file
	while(my $seq_obj = $input_fastq->next_seq()) {
		my $original_fasta = $seq_obj->seq();
		
		# lets get a random extraction so many times as DNA copies were chosen
		for my $i (1 .. $$opts{dna_copies}) {
			
			if ($$opts{verbose}) {print "Extracting DNA copy $i\n";}
			
			my $seq_fasta = $original_fasta;
			my $dna_size = 0;
			# As the sequence is trimmed the fragments are saved in an array
			while(length($seq_fasta) >= $dna_size) {
				# fragments size are obtained according to a Gaussian distribution
				$dna_size = int(random_normal(1,$$opts{extraction_size},$$opts{std_dev}));
				
				# print "dna size: $dna_size, fasta length: ".length($seq_fasta)."\n";
				
				push(@dna,substr($seq_fasta,0,$dna_size));
				if ($dna_size < length($seq_fasta)) {
					$seq_fasta = substr($seq_fasta,$dna_size,length($seq_fasta));
				}
			}
		}
	}
	
	if ($$opts{verbose}) {print "\nfragments number: ".(scalar($#dna)+1)."\n";}
	
	return(\@dna);
}

# Lets simulate the distribution of DNA through the aliquots
sub dna_to_aliquots {
	my $opts = shift;
	my $dna = shift;
	my %plate96;
	
	# read the tags file for saving them in a Hash of arrays
	open (my $aliq_fh, $$opts{aliquot_tags}) || die ("\nERROR: the file ".$$opts{aliquot_tags}." could not be found\n");
	while (my $line = <$aliq_fh>) {
		chomp($line);
		$plate96{$line}=();
	}
	close ($aliq_fh);
	
	if ($$opts{verbose}) {print "Tags loaded from file ".$$opts{aliquot_tags}."\n";}
	
	# total number of fragments of DNA in pool
	my $fragments_num = int((scalar($#{$dna})+1)/$$opts{dna_copies});
	# number of fragments to put in each aliquot according with the chosen genome percentage
	my $amount = int($$opts{genome_percent}*$fragments_num);
	my $random_index = 0;
	
	if ($$opts{verbose}) {
		print "fragments number to get the genome size: $fragments_num, ".$$opts{genome_percent}."% --> $amount\n";
	}
	
	# Lets load random DNA sequences in the array of each aliquot until the amount fixed
	foreach my $tag (keys %plate96) {
		for my $j (0..$amount-1) {
			$random_index = int(rand(scalar($#{$dna})));
			push(@{$plate96{$tag}},$$dna[$random_index]);
		}
	}
	@$dna=();
	return(\%plate96);
}

# get total digested DNA
sub total_digestion{
	my $opts = shift;
	my $aliq_HoA = shift;
	
	my @dna;
	my %digested_aliq;
	
	my $rs = $$opts{total_rs};
	$rs =~ s/\^//;

	foreach my $tag (keys %$aliq_HoA) {
		for my $j (0..$#{$$aliq_HoA{$tag}}) {
			
			my $seq_fasta = $$aliq_HoA{$tag}[$j];
			my @fragments = split($rs,$seq_fasta);
		
			push(@{$digested_aliq{$tag}},@fragments);
		}
	}
	%$aliq_HoA=();
	
	if ($$opts{verbose}) {
		my $sum_of_rs = 0;
		foreach my $tag (keys %digested_aliq) {
			$sum_of_rs += (scalar($#{$digested_aliq{$tag}})+1);
		}
		print "Restriction sites found on average per aliquot: ".int($sum_of_rs/scalar(keys %digested_aliq))."\n";
	}
	return(\%digested_aliq);
}

sub print_happy_reads {
	my $opts = shift;
	my $aliq_HoA = shift;
	
	my ($rs1,$rs2) = split(/\^/,$$opts{total_rs});
	my $read_length = ($$opts{reads_length} - length($rs2));
	
	# Lets print the reads
	open (my $o_fh, ">output/happy_reads.fasta") || die ("\nERROR: the file 'output/happy_reads.fasta' could not be found\n");
	
	my $counter = 0;
	foreach my $tag (keys %$aliq_HoA) {
		$counter = 0;
		for my $j (0..$#{$$aliq_HoA{$tag}}) {
			my $read = $$aliq_HoA{$tag}[$j];
			
			if (length($read) >= ($$opts{reads_length} - length($rs2))) {
				$counter++;
				my $trimmed_read = substr($read,0,$read_length);
				print $o_fh ">${tag}_read${counter}\n${rs2}${trimmed_read}\n";
			}
		}
	}
}

################################################################################
#                                    MAIN                                      #
################################################################################
if ($Options{verbose}) {print "---------------Extracting DNA\n";}
# ------------------------------------------------------------------------------ DNA extraction
my $dna = dna_extraction(\%Options);

if ($Options{verbose}) {print "---------------Distributing DNA in aliquots\n";}
# ------------------------------------------------------------------------------ Distributing a percentage of genome in aliquots
my $aliquots = dna_to_aliquots(\%Options, $dna);

if ($Options{verbose}) {print "---------------Digesting DNA\n";}
# ------------------------------------------------------------------------------ total digestion
my $digested_dna = total_digestion(\%Options, $aliquots);

if ($Options{verbose}) {print "---------------printing happy reads\n";}
# ------------------------------------------------------------------------------ print results
print_happy_reads(\%Options, $digested_dna);






