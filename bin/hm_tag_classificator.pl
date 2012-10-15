#!/usr/bin/perl

=head1 NAME

 hm_tag_classificator
 This script classifies markers from happy mapping by adding the aliquot tag to their name.
 sequence tags are trimmed until the restriction site and 
 markers without tag or restriction site are discarded.
 Output reads are trimmed to a minimum fixed length.

=cut

=head1 AUTHORS

  Noe Fernandez-Pozo
  (nf232@cornell.edu)

=cut

use strict;
use warnings;

use Getopt::Long;
use Bio::SeqIO;
use threads;

################################################################################
#                               LOAD & CHECKING                                #
################################################################################


if (!$ARGV[0]) {
	help_msg();
	exit;
}

# 12 Oct 2012 -> kepp RS in markers seq
my $VERSION = "0.0.2";

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
	$h{seqs} = 10000;
	$h{cores} = 1;
	$h{r_site} = 'TGCAG';
	
	GetOptions( \%h,
		"seqs|s=i",
		"cores|c=i",
		"fastq|f=s",
		"tags|t=s",
		"help|h",
		"length|l=i",
		"r_site|r=s",
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
	if ((!$options{fastq}) || (!$options{tags})) {
		print "\nERROR: a fasta file and a list of tags must be provided\n\n";
		
		help_msg();
		exit;
	}
	
	# -------------------------------------------- check files
	check_file_exist($options{fastq},"File ".$options{fastq}." was not found. A file with the marker reads is required");
	check_file_exist($options{tags},"File ".$options{tags}." was not found. A file with tags used for each aliquot of the happy mapping is required:\n\nacag\nactg\nagac\n...\n");
	
	if ($options{verbose}) {
		print STDERR "--- YOUR PARAMETERS AND OPTIONS WERE CORRECTLY CHECKED ---\n";
	}
	
	# -------------------------------------------- removing old files and making new directories
	if (-e "tmp_files/"){
		`rm -rf tmp_files`;
	}
	if (-e "output/") {
		`rm -rf output`;
	}
	mkdir ("output");
	mkdir ("tmp_files");
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

Usage: $0 -f <file.fastq> -t <tags_file.txt> [other options or parameters]

  *** MANDATORY OPTIONS ****************************************************************
  
  -fastq|f <file>           A fastq file with the happy mapping reads
  -tags|t <file>            A text file with the happy mapping tags
 
  *** MORE OPTIONS AND PARAMETERS ******************************************************

  -length|l <integer>       Minimum length of the output markers (default=69bp)
  -r_site|r <string>        Sequence expected for the restriction site (default='TGCAG')
  -cores|c <integer>        Number of cores used in the analysis (default=1)
  -seqs|s <integer>         Number sequences sent to each core (default=10000)
  
  -verbose                  Verbose mode
  -version|v                Shows program version
  -help|h                   Shows this help

END

}


################################################################################
#                                 SUBROUTINES                                  #
################################################################################

# loading tags from Happy Tags file in a hash
sub load_tags{
	
	my $tags_file = shift;
	my %tags;
	
	open (my $tag_fh, $tags_file) || die ("\nERROR: the file '$tags_file' could not be found\n");
	while (my $line = <$tag_fh>) {
		chomp($line);
		$tags{uc($line)} = 1;
	}
	return \%tags;
}

# split fasta in chunks and save it into an array for parrallel execution
sub split_fasta {
	
	my $opt_h = shift;
	my $tag_h = shift;
	
	my %opts = %$opt_h;
	my %tags = %$tag_h;
	
	my $threads = $opts{cores};
	my $min_length = $opts{length};
	my $r_site = $opts{r_site};
	my $fasta_file = $opts{fastq};
	my $v = $opts{verbose};
	
	my $chunk_lines = ($opts{seqs}*4);
	my $line_counter = 0; # useful to create the chunks
	my $counter = 0; # to create fasta tmp files
	my @chunk; # for saving the chunks
	my @array_thr; #for saving the threads
	my $string_chunk = '';
	my @thr_args;
	my $running_thrs = 0;
	
	open (my $thread_fh, ">tmp_files/tmp_threads_file.txt") || die ("ERROR: was not possible to create 'tmp_files/tmp_threads_file.txt'");
	open (my $error_fh, ">output/seqs_without_tag.txt") || die ("ERROR: was not possible to create 'output/seqs_without_tag.txt'");
	open (my $fq_fh, $fasta_file) || die ("\nERROR: the file '$fasta_file' could not be found\n");
	
	while (my $line = <$fq_fh>) {
		chomp($line);
		$line_counter++;
		push(@chunk,$line);
		
		# when a chunk is filled
		if ($line_counter % $chunk_lines == 0){
			close_thread(\@array_thr);
			$counter++;
			
			my @running = threads->list(threads::running);
			$running_thrs = (scalar($#running)+1);
			
			$string_chunk = join("\n",@chunk);
			@thr_args=($string_chunk,\%tags,$r_site,$counter,$min_length);
			
			# do not start more threads until we finished the ones we are using
			while ($running_thrs >= $threads) {
				sleep(1);
				@running = threads->list(threads::running);
				$running_thrs = (scalar($#running)+1);
			}
			
			my $thr = threads->create('one_thread', \@thr_args);
			push(@array_thr, $thr);
			
			if($v) {
				print "threads running: ".($running_thrs+1).", seqs loaded: ".($line_counter/4)."\n";
			}
			
			@chunk = ();
			$string_chunk = '';
		}
	}
	# ---------------------------------------------rest of the sequences are analised
	$string_chunk = join("\n",@chunk);
	$counter++;
	@thr_args=($string_chunk,\%tags,$r_site,$counter,$min_length);
	my $thr = threads->create('one_thread', \@thr_args);
	push(@array_thr, $thr);
	close_thread(\@array_thr);
	
	if($v) {
		print "threads running: $running_thrs, seqs loaded: ".($line_counter/4)."\n";
	}
	# -------------------------------------------------------------------------------
	
	# wait until all the threads are finished
	foreach my $thr (@array_thr) {
		while ($thr->is_running()) {
			sleep(1);
		}
		if ($thr->is_joinable()) {
			$thr->join();
		}
	}
	
	# close files
	close($thread_fh);
	close($error_fh);
	
	return ($line_counter);
}

# analysis of sequences in each thread
sub one_thread {
	my $input_argv = shift;
	my @argv_a = @$input_argv;
	
	my ($string_fasta, $tags_h, $r_site, $pid, $min_length) = @argv_a;
	
	my %tags = %$tags_h;
	my $tag_out_fh='';
	my %saved_seqs;
	my %no_tag_seqs;
	my %tag_counter;
	my $no_tag_counter = 0;
	my $too_short_counter = 0;
	
	# get sequences one by one and analyse them
	open(my $seqs_ah, "<", \$string_fasta);
	my $input_fastq = Bio::SeqIO->new(-fh => $seqs_ah, -format => 'fastq');
	$tag_out_fh = Bio::SeqIO->new( -file => ">tmp_files/preprocessed_$pid.fasta", -format => 'fasta');
	
	while(my $seq_obj = $input_fastq->next_seq()) {
		
		my $seq_name = $seq_obj->id();
		my $seq_fasta = $seq_obj->seq();
		
		# sequence with a tag before the restriction site
		if ($seq_fasta =~ m/([a-z]{4})${r_site}(.+)/i){
			my $tag = $1;
			my $rest_of_Seq = $2;
			
			# sequences with a right tag (contained in our happy tags file)
			if ($tags{uc($tag)}) {
				$seq_name = "${tag}_${seq_name}";
				if ($seq_fasta =~ m/(.*${tag}${r_site})/i){
					$seq_fasta =~ s/.*${tag}${r_site}/$r_site/i;
				}
				if (length($seq_fasta) >= $min_length) {
					my $seq = substr($seq_fasta,0,$min_length);
					my $new_seq = Bio::Seq->new(-id => $seq_name, -seq => $seq);
					
					$tag_counter{$tag}++;
					$tag_out_fh->write_seq($new_seq);
					$saved_seqs{$seq_name} = 1;
					
				# sequences with tag but too short
				} else {
					$no_tag_seqs{$seq_name} = "shorter than $min_length bp";
					$too_short_counter++;
				}
			} else { # first tag found was not right
				# there is a second tag
				if ($rest_of_Seq =~ /([a-z]{4})${r_site}/i){
					my $tag2 = $1;
					
					if ($tags{uc($tag2)}) { # second tag is right
						$seq_name = "${tag2}_${seq_name}";
						if ($seq_fasta =~ m/(.*${tag2}${r_site})/i){
							$seq_fasta =~ s/.*${tag2}${r_site}/$r_site/i;
						}
						
						if (length($seq_fasta) >= $min_length) {
							my $seq = substr($seq_fasta,0,$min_length);
							my $new_seq = Bio::Seq->new(-id => $seq_name, -seq => $seq);

							$tag_counter{$tag2}++;
							$tag_out_fh->write_seq($new_seq);
							$saved_seqs{$seq_name} = 1;
							
						} else {
							$no_tag_seqs{$seq_name} = "shorter than $min_length bp in 2nd tag";
							$too_short_counter++;
						}
					} else { # 2nd tag found was not right
						$no_tag_seqs{$seq_name} = "2nd tag found was wrong: $tag";
						$no_tag_counter++;
					}
				} else {
					$no_tag_seqs{$seq_name} = "the tag found was wrong: $tag";
					$no_tag_counter++;
				}
			}
			
		# sequences without tag
		} else {
			$no_tag_seqs{$seq_name} = "no tag found";
			$no_tag_counter++;
		}
	}
	open (my $thread_fh, ">>tmp_files/tmp_threads_file.txt");
	open (my $error_fh, ">>output/seqs_without_tag.txt");
	
	print $thread_fh "no_tag_num: $no_tag_counter\n";
	print $thread_fh "too_short_num: $too_short_counter\n";
	print $thread_fh "with_tag_num: ".(keys %saved_seqs)."\n";
	
	foreach my $tag (keys %tag_counter){
		print $thread_fh "$tag: $tag_counter{$tag}\n";
	}
	foreach my $seq_name (keys %no_tag_seqs){
		print $error_fh "$seq_name\t$no_tag_seqs{$seq_name}\n";
	}
	
}

# close threads with finished jobs
sub close_thread {
	my $thr_a = shift;
	my @array_thr = @$thr_a;
	
	foreach my $thr (@array_thr) {
		if ($thr->is_joinable()) {
			$thr->join();
		}
	}
}

sub print_stats {
	
	my $fasta_lines = shift;
	my $verbose = shift;
	
	my $total_seqs = ($fasta_lines/4);
	my $total_tagged = 0;
	my $total_no_tagged = 0;
	my $total_too_short = 0;
	my %total_tag_counter;

	open (my $result_fh, "tmp_files/tmp_threads_file.txt") || die ("\nERROR: the file 'tmp_files/tmp_threads_file.txt' could not be found\n");
	while (my $line = <$result_fh>) {
	
		if ($line =~ /^no_tag_num:\s(\d+)$/) {
			$total_no_tagged += $1;
		}
		elsif ($line =~ /^too_short_num:\s(\d+)$/) {
			$total_too_short += $1;
		}
		elsif ($line =~ /^with_tag_num:\s(\d+)$/){
			$total_tagged += $1;
		}
		elsif ($line =~ /^([acgt]+):\s(\d+)$/i){
			$total_tag_counter{$1} += $2;
		}
	}
	
	# Lets print some stats
	open (my $o_fh, ">output/stats.txt") || die ("\nERROR: the file 'output/stats.txt' could not be found\n");
	
	
	# print aliquots stats
	print $o_fh "\nNumber of markers found by aliquot:\n";
	foreach my $tag (keys %total_tag_counter){
		print $o_fh "$tag: $total_tag_counter{$tag}\n";
	}

	print $o_fh "\nNumber of sequences: $total_seqs\n";
	print $o_fh "total tagged: $total_tagged (".sprintf("%.2f", (100*$total_tagged/$total_seqs))."%)\n";
	print $o_fh "total NO tagged: $total_no_tagged (".sprintf("%.2f", (100*$total_no_tagged/$total_seqs))."%)\n";
	print $o_fh "total too short: $total_too_short (".sprintf("%.2f", (100*$total_too_short/$total_seqs))."%)\n";

	if ($verbose) {
		print "\nNumber of sequences: $total_seqs\n";
		print "total tagged: $total_tagged (".sprintf("%.2f", (100*$total_tagged/$total_seqs))."%)\n";
		print "total NO tagged: $total_no_tagged (".sprintf("%.2f", (100*$total_no_tagged/$total_seqs))."%)\n";
		print "total too short: $total_too_short (".sprintf("%.2f", (100*$total_too_short/$total_seqs))."%)\n";
	
		print "\ntotal analysed: $total_tagged + $total_no_tagged + $total_too_short \=".($total_tagged+$total_no_tagged+$total_too_short)."\n\n";
	}
}


################################################################################
#                                    MAIN                                      #
################################################################################

# ------------------------------------------------------------------------------ loading tags
my $tags_hash = load_tags($Options{tags});
my %tags = %$tags_hash;

# ------------------------------------------------------------------------------ split fasta and send chunks to several cores
my ($line_counter) = split_fasta(\%Options, \%tags);

# ------------------------------------------------------------------------------ reading and printing stats
print_stats($line_counter, $Options{verbose});

# ------------------------------------------------------------------------------ join temporal fasta files and remove temporal file
`cat tmp_files/preprocessed_*.fasta > output/labeled_seqs.fasta`;
`rm -rf tmp_files`;

