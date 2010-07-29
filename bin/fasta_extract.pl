#!/usr/bin/perl

=head1 NAME

 fasta_extract.pl
 Tool to extract sequences from a fasta file in different ways

=cut

=head1 SYPNOSIS

 fasta_extract.pl [-h] -f <input_fasta_file> [-o output_basename]
                       [-i <extract_by_id>] [-l <extract_by_length>] [-r <extract_by_regexp>] 
                       [-v]

=head2 I<Flags:>

=over


=item -f

B<input_fasta_file>       input fasta file (mandatory)

=item -o

B<output_basename>        output basename (by default it will printed as stdout)

=item -i

B<extract_by_id>          filename to extract sequences by id

=item -l 

B<extract_by_length>      integer to extract sequences by length     

=item -r

B<regexp>                 regular expression to extract the sequences

=item -z

B<format>                 format

=item -v

B<invert_selection>       invert the selection

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This script extract a list of sequences from a fasta/qual file. 

 There are three ways to supply the list of sequences:

   1) As list of ids in a file (-i <id_list_file>)

   2) As length value (-l <length>)

   3) As regular expression to match with the id (-r <regexp>)

   Examples:

      fasta_extract.pl -f test.fasta -i id_file.txt > extract.fasta

      fasta_extract.pl -f test.fasta -l 100 -v > seq_less100.fasta

      fasta_extract.pl -f test.fasta -r 'Sl' > sl_seq.fasta

      fasta_extract.pl -f test.fasta -r 'AAAAAAAAAA$' > seq_with_poliA.fasta

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 fasta_extract.pl


=cut

use strict;
use warnings;

use Getopt::Std;
use Bio::SeqIO;

our ($opt_f, $opt_o, $opt_i, $opt_l, $opt_r, $opt_v, $opt_z, $opt_h);
getopts("f:o:i:l:r:vz:h");
if (!$opt_f && !$opt_o && !$opt_i && !$opt_l && !$opt_r && !$opt_v && !$opt_z && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}

## Get the arguments and check them

my $infasta_file = $opt_f || 
    die("INPUT ARGUMENTS ERROR: Input option was not supplied (-f <input_fasta_file>).\n");

my $id_file = $opt_i;
my $length_cf = $opt_l;
my $regexp;
if (defined $opt_r) {
    $regexp = qr/$opt_r/;
}

my $format = $opt_z || 'fasta';

if (!$id_file && !$length_cf && !$regexp) {
    die("EXTRACTION ARGUMENTS ERROR: None of the extract options were used (-i <extract_by_id>,-l <extract_by_length>,-r <regexp>)\n");
}

my $out = $opt_o;

## Extract the list of ids from the id file if exists

my %ids;

if (defined $id_file) {

    print STDERR "\nEXTRACT BY IDS option [Enabled].\n";

    open my $idfh, '<', $id_file || 
	die("OPEN FILE ERROR: File=$id_file can not be openned (system error:$!)\n");

    while (<$idfh>) {
	chomp($_);

	## It will take only the first column and remove the rest

	my @data = split(/\t/, $_);

	## It will ignore the duplications overwritting the ids

	$ids{$data[0]} =1;
    }

    my $ids_count = scalar(keys %ids);
    print STDERR "\t$ids_count ids were parsed from $id_file.\n"
}

## Open the input and the output fasta file


print STDERR "\nEXTRACTING SEQUENCES FROM $infasta_file.\n";

my $c = 0;

my $inio = Bio::SeqIO->new(
                             -file   => $infasta_file,
                             -format => $format
                          );


## It will print the sequences in STDOUT by default

my $outio;
if (defined $out) {

    my $outfasta = $out . '.fasta';

    $outio = Bio::SeqIO->new( 
	                      -file   => ">$outfasta",
	                      -format => $format
	                    );	
}
else {
    $outio = Bio::SeqIO->new( -format => $format );
}

my $s = 0;
while (my $seqobj = $inio->next_seq()) {
    
    $s++;
    my $id = $seqobj->primary_id();
    my $lg = $seqobj->length();

    print STDERR "Processing id=$id (sequence $s)    \r";

    ## Select only if it have one of the extract option

    my $sel = 0;
	
    if (defined $id_file) {
	if (exists $ids{$id}) {
	    $sel = 1;
	}
    }
    if (defined $length_cf) {
	if ($lg > $length_cf) {
	    $sel = 1;
	}
    }
    if (defined $regexp) {
	if ($id =~ $regexp) {
	    $sel = 1;
	}
	elsif ($seqobj->seq() =~ $regexp) {
	    $sel = 1;
	}
    }

    if ($opt_v) {
	
	## Invert the selection 
	$sel =~ s/0/2/;
	$sel =~ s/1/0/;
	$sel =~ s/2/1/;
    }

    ## If the sequence was selected, it will printed

    if ($sel == 1) {
	$c++;
	if ($format eq 'fastq') {
	    $outio->write_fastq($seqobj);
	}
	else {
	    $outio->write_seq($seqobj);
	}
    }
}

print STDERR "\n\nDONE.\n\t$c sequences were extracted from the $infasta_file.\n\n";



=head2 help

  Usage: help()
  Desc: print help of this script
  Ret: none
  Args: none
  Side_Effects: exit of the script
  Example: if (!@ARGV) {
               help();
           }

=cut

sub help {
  print STDERR <<EOF;
  $0:

    Description:

      This script extract a list of sequences from a fasta/qual file. 

      There are three ways to supply the list of sequences:

	  1) As list of ids in a file (-i <id_list_file>)
    
          2) As length value (-l <length>)

          3) As regular expression to match with the id (-r <regexp>)

    Usage:
     
      fasta_extract.pl [-h] -f <input_fasta_file> [-o <output_basename>] [-i <filter_by_id>] [-l <filter_by_length>] 
	                                          [-r <filter_by_regexp>] [-v]
      
    Examples:

      fasta_extract.pl -f test.fasta -i id_file.txt > extract.fasta

      fasta_extract.pl -f test.fasta -l 100 -v > seq_less100.fasta

      fasta_extract.pl -f test.fasta -r 'Sl' > sl_seq.fasta

      fasta_extract.pl -f test.fasta -r 'AAAAAAAAAA$' > seq_with_poliA.fasta

    Flags:

      -f <input_fasta_file>       input fasta file (mandatory)
      -o <output_basename>        output basename (by default it will printed as stdout)
      -i <extract_by_id>          filename to extract sequences by id
      -l <extract_by_length>      integer to extract sequences by length     
      -r <regexp>                 regular expression to extract the sequences
      -v <invert_selection>       invert the selection
      -z <format>                 input file format  
      -h <help>                   print the help
     

EOF
exit (1);
}

