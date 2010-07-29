#!/usr/bin/perl

=head1 NAME

 bio_fileformat_tool.pl
 This script change the format betwen files using the Bio::SeqIO (version.1.0.).

=cut

=head1 SYPNOSIS

 bio_fileformat_tool.pl [-h] -i <input_file> -a <input_format> -o <output_file> -z <output_format> -f <filter_file> -R

=head2 I<Flags:>

=over


=item -i

B<input_file>             input file (mandatory)

=item -a

B<input_format>           input format (mandatory)

=item -o

B<output_file>            output file (mandatory)

=item -z

B<output_format>          output format (mandatory)

=item -f 

B<filter_file>            a list of id to exclude in the new file (it will ignore everything after the first column)

=item -R

B<reverse>                use a reverse option with the filter and select only the ids of the filter file

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This script change the format of the sequence file using the Bio::SeqIO options. The format supported are:

  Also can be used to filter (-f) or extract (-f -R) a list of ids from a file. 

+==================+==========================================+=================================+
| Name             | Description                              | File extension                  |
+==================+==========================================+=================================+
| abi              | ABI tracefile                            | ab[i1]                          |
+------------------+------------------------------------------+---------------------------------+
| ace              | Ace database                             | ace                             |
+------------------+------------------------------------------+---------------------------------+
| agave            | AGAVE                                    | XML                             |  
+------------------+------------------------------------------+---------------------------------+ 	
| alf  	           | ALF tracefile                            | alf  	                        |
+------------------+------------------------------------------+---------------------------------+
| asciitree        | write-only, to visualize features        |                                 |
+------------------+------------------------------------------+---------------------------------+
| bsml 	           | BSML, using XML::DOM                     | bsml                            |
+------------------+------------------------------------------+---------------------------------+
| bsml_sax         | BSML, using XML::SAX                     |		                        |
+------------------+------------------------------------------+---------------------------------+
| chadoxml         | CHADO sequence format                    | 	                        |
+------------------+------------------------------------------+---------------------------------+
| chaos            | CHAOS sequence format                    | 		                |
+------------------+------------------------------------------+---------------------------------+
| chaosxml         | Chaos XML                                |		                        |
+------------------+------------------------------------------+---------------------------------+
| ctf 	           | CTF tracefile 	                      | ctf 	                        |  
+------------------+------------------------------------------+---------------------------------+
| embl 	           | EMBL database 	                      | ebl | emb | dat                 |
+------------------+------------------------------------------+---------------------------------+
| entrezgene       | Entrez Gene                              | ASN1   		                |
+------------------+------------------------------------------+---------------------------------+
| excel            | Excel 		                      |                                 |
+------------------+------------------------------------------+---------------------------------+
| exp 	           | Staden EXP format                        | exp    	                        |
+------------------+------------------------------------------+---------------------------------+
| fasta            | FASTA 	                              | fast | seq | fa | fsa | nt | aa |
+------------------+------------------------------------------+---------------------------------+
| fastq            | quality score data in FASTA-like format  | fastq                           |
+------------------+------------------------------------------+---------------------------------+
| flybase_chadoxml | variant of Chado XML                     |                                 |
+------------------+------------------------------------------+---------------------------------+ 	   
| game 	           | GAME XML 	                              |                                 |
+------------------+------------------------------------------+---------------------------------+
| gcg 	           | GCG 	                              | gcg 	                        |
+------------------+------------------------------------------+---------------------------------+
| genbank 	   | GenBank 	                              | gbank | genbank                 |
+------------------+------------------------------------------+---------------------------------+
| interpro 	   | InterProScan XML                         | 	                        |
+------------------+------------------------------------------+---------------------------------+
| kegg 	           | KEGG 		                      |                                 |
+------------------+------------------------------------------+---------------------------------+
| largefasta 	   | Large files, fasta format                | 		                |
+------------------+------------------------------------------+---------------------------------+
| lasergene 	   | Lasergene format                         | 		                |
+------------------+------------------------------------------+---------------------------------+
| locuslink 	   | LocusLink LL_tmpl                        | 	                        |
+------------------+------------------------------------------+---------------------------------+
| metafasta 	   |                                          |                                 |
+------------------+------------------------------------------+---------------------------------+
| phd 	           | Phred 	                              | phred                           |
+------------------+------------------------------------------+---------------------------------+
| pir 	           | PIR database 	                      | pir 	                        |
+------------------+------------------------------------------+---------------------------------+
| pln 	           | PLN tracefile                            | pln                             |
+------------------+------------------------------------------+---------------------------------+
| qual             | Phred (quality scores)                   | 	                        |
+------------------+------------------------------------------+---------------------------------+
| raw 	           | plain text 	                      | txt 	                        |
+------------------+------------------------------------------+---------------------------------+
| scf 	           | Standard Chromatogram Format             | scf                             |
+------------------+------------------------------------------+---------------------------------+
| strider          | DNA Strider format                       | 	                        |
+------------------+------------------------------------------+---------------------------------+
| swiss 	   | SwissProt 	                              | sp 	                        |
+------------------+------------------------------------------+---------------------------------+
| tab 	           | tab-delimited                            |  		                |
+------------------+------------------------------------------+---------------------------------+
| table 	   | Table 		                      |                                 |
+------------------+------------------------------------------+---------------------------------+
| tigr 	           | TIGR XML                                 | 	                        |
+------------------+------------------------------------------+---------------------------------+
| tigrxml 	   | TIGR Coordset XML                        | 	                        |
+------------------+------------------------------------------+---------------------------------+
| tinyseq 	   | NCBI TinySeq XML                         | 		                |
+------------------+------------------------------------------+---------------------------------+
| ztr 	           | ZTR tracefile                            | ztr                             |
+==================+==========================================+=================================+
| blastm8          | Blast results in m8 format               | m8                              |
+==================+==========================================+=================================+


=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 bio_fileformat_tool.pl


=cut

use strict;
use warnings;

use Getopt::Std;
use Bio::SeqIO;

our ($opt_i, $opt_a, $opt_o, $opt_z, $opt_f, $opt_R, $opt_h);
getopts("i:a:o:z:f:Rh");
if (!$opt_i && !$opt_a && !$opt_o && !$opt_z && !$opt_f && !$opt_R && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}

my $input_file = $opt_i || die("DATA ARGUMENT -i input_file WAS NOT SUPPLIED.\n");
my $input_format = $opt_a || die("DATA ARGUMENT -a input_format WAS NOT SUPPLIED.\n");
my $output_file = $opt_o || die("DATA ARGUMENT -o output_file WAS NOT SUPPLIED.\n");
my $output_format = $opt_z || die("DATA ARGUMENT -z output_format WAS NOT SUPPLIED.\n");

## Use the blastm8 restriction, only to filter this blast files

if ($input_format eq 'blastm8' && $output_format ne 'blastm8') {
    die("DATA ARGUMENT ERROR: Blastm8 format only can be used to be filtered, it can not be changed to other formats.\n");
}
elsif ($input_format ne 'blastm8' && $output_format eq 'blastm8') {
    die("DATA ARGUMENT ERROR: Blastm8 format only can be used to be filtered, it can not be changed to other formats.\n");
}





my %excl_list;

if ($opt_f) {

    print STDERR "\nPARSING THE FILTER FILE:\n\n";

    open my $filter_fh, '<', $opt_f ||
	die("Sorry, I can not open the filter file: $opt_f.\n");
    
    my $t = 0; 
    while (<$filter_fh>) {
	chomp($_);
	$t++;
	
	## This only will take the first column
	my @data = split(/\t/, $_);

	print STDERR "\tProcessing line: $t with id=$data[0] for filter file                      \r";

	## Just ignore the duplications

	$excl_list{$data[0]} = 1;
    }
    print STDERR "\n\tFilter file processed.\n";

}


my $s = 0;
my $n = 1;

if ($input_format ne 'blastm8') {
   
    my $out = Bio::SeqIO->new( 
                               -file => ">$output_file",
                               -format => $output_format
                             );

    print STDERR "\nPROCESSING FILE:$input_file FROM FORMAT:$input_format TO:$output_format.\n\n";

    
    my $in  = Bio::SeqIO->new(
	                       -file   => $input_file,
                               -format => $input_format
                             );
 
    while ( my $seq = $in->next_seq() ) {
	my $id = $seq->display_id();
	    
	print STDERR "\tProcesing $n sequence/qscore/id ($id)                 \r";
	
	if ($opt_f) {
	    
	    ## For filter purposes, we remove some characters from the id
	    $id =~ s/;//;
	    
	    if ($opt_R) {
		if (exists $excl_list{$id}) {
		    $out->write_seq($seq);
		    $s++;
		}
	    }
	    else {
		unless (exists $excl_list{$id}) {
		    $out->write_seq($seq);
		    $s++;
		}
	    }
	}
	else {
	    $out->write_seq($seq);
	    $s++;
	} 
	$n++;
    }
}
else {

    open my $in, '<', $input_file ||
	die ("Sorry, I can not open the input file $input_file.\n");

    open my $out, '>', $output_file ||
	die ("Sorry, I can not open the output file $output_file.\n");

    while(<$in>) {
	chomp($_);
	my @indata = split(/\t/, $_);
	my $id = $indata[0];

	print STDERR "\tProcesing $n sequence/qscore/id ($id)                 \r";
	
	if ($opt_f) {
	    
	    ## For filter purposes, we remove some characters from the id
	    $id =~ s/;//;
	    
	    if ($opt_R) {
		if (exists $excl_list{$id}) {
		    print $out "$_\n";
		    $s++;
		}
	    }
	    else {
		unless (exists $excl_list{$id}) {
		    print $out "$_\n";
		    $s++;
		}
	    }
	}
	else {
	    print $out "$_\n";
	    $s++;
	} 
	$n++;

    }
}
print STDERR "\n\nDONE, $s seq/qscores/ids were copied to the file:$output_file in $output_format format.\n\n";


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
       This script change the format of the sequence file using the Bio::SeqIO options. The format supported are:

         +==================+==========================================+=================================+
         | Name             | Description                              | File extension                  |
         +==================+==========================================+=================================+ 
         | abi              | ABI tracefile                            | ab[i1]                          |
         +------------------+------------------------------------------+---------------------------------+
         | ace              | Ace database                             | ace                             |
         +------------------+------------------------------------------+---------------------------------+
         | agave            | AGAVE                                    | XML                             |  
         +------------------+------------------------------------------+---------------------------------+ 	
         | alf              | ALF tracefile                            | alf  	                         |
         +------------------+------------------------------------------+---------------------------------+
         | asciitree        | write-only, to visualize features        |                                 |
         +------------------+------------------------------------------+---------------------------------+
         | bsml             | BSML, using XML::DOM                     | bsml                            |
         +------------------+------------------------------------------+---------------------------------+
         | bsml_sax         | BSML, using XML::SAX                     |	                         |
         +------------------+------------------------------------------+---------------------------------+
         | chadoxml         | CHADO sequence format                    | 	                         |
         +------------------+------------------------------------------+---------------------------------+
         | chaos            | CHAOS sequence format                    | 		                 |
         +------------------+------------------------------------------+---------------------------------+
         | chaosxml         | Chaos XML                                |	                         |
         +------------------+------------------------------------------+---------------------------------+
         | ctf 	            | CTF tracefile 	                       | ctf 	                         |  
         +------------------+------------------------------------------+---------------------------------+
         | embl             | EMBL database 	                       | ebl | emb | dat                 |
         +------------------+------------------------------------------+---------------------------------+
         | entrezgene       | Entrez Gene                              | ASN1   	                 |
         +------------------+------------------------------------------+---------------------------------+
         | excel            | Excel 		                       |                                 |
         +------------------+------------------------------------------+---------------------------------+
         | exp 	            | Staden EXP format                        | exp    	                 |
         +------------------+------------------------------------------+---------------------------------+
         | fasta            | FASTA 	                               | fast | seq | fa | fsa | nt | aa |
         +------------------+------------------------------------------+---------------------------------+
         | fastq            | quality score data in FASTA-like format  | fastq                           |
         +------------------+------------------------------------------+---------------------------------+
         | flybase_chadoxml | variant of Chado XML                     |                                 |
         +------------------+------------------------------------------+---------------------------------+ 	   
         | game             | GAME XML 	                               |                                 |
         +------------------+------------------------------------------+---------------------------------+
         | gcg 	            | GCG 	                               | gcg 	                         |
         +------------------+------------------------------------------+---------------------------------+
         | genbank 	    | GenBank 	                               | gbank | genbank                 |
         +------------------+------------------------------------------+---------------------------------+
         | interpro 	    | InterProScan XML                         | 	                         |
         +------------------+------------------------------------------+---------------------------------+
         | kegg             | KEGG 		                       |                                 |
         +------------------+------------------------------------------+---------------------------------+
         | largefasta 	    | Large files, fasta format                | 		                 |
         +------------------+------------------------------------------+---------------------------------+
         | lasergene 	    | Lasergene format                         | 		                 |
         +------------------+------------------------------------------+---------------------------------+ 
         | locuslink 	    | LocusLink LL_tmpl                        | 	                         |
         +------------------+------------------------------------------+---------------------------------+
         | metafasta 	    |                                          |                                 |
         +------------------+------------------------------------------+---------------------------------+
         | phd 	            | Phred 	                               | phred                           |
         +------------------+------------------------------------------+---------------------------------+
         | pir 	            | PIR database     	                       | pir 	                         |
         +------------------+------------------------------------------+---------------------------------+
         | pln 	            | PLN tracefile                            | pln                             |
         +------------------+------------------------------------------+---------------------------------+
         | qual             | Phred (quality scores)                   | 	                         |
         +------------------+------------------------------------------+---------------------------------+
         | raw 	            | plain text 	                       | txt 	                         |
         +------------------+------------------------------------------+---------------------------------+
         | scf 	            | Standard Chromatogram Format             | scf                             |
         +------------------+------------------------------------------+---------------------------------+
         | strider          | DNA Strider format                       | 	                         |
         +------------------+------------------------------------------+---------------------------------+
         | swiss 	    | SwissProt 	                       | sp 	                         |
         +------------------+------------------------------------------+---------------------------------+
         | tab 	            | tab-delimited                            |  		                 |
         +------------------+------------------------------------------+---------------------------------+
         | table 	    | Table 		                       |                                 |
         +------------------+------------------------------------------+---------------------------------+
         | tigr             | TIGR XML                                 | 	                         |
         +------------------+------------------------------------------+---------------------------------+
         | tigrxml 	    | TIGR Coordset XML                        | 	                         |
         +------------------+------------------------------------------+---------------------------------+
         | tinyseq 	    | NCBI TinySeq XML                         | 		                 | 
         +------------------+------------------------------------------+---------------------------------+
         | ztr 	            | ZTR tracefile                            | ztr                             |
         +==================+==========================================+=================================+
	 | blastm8          | Blast results in m8 format               | m8                              |
	 +==================+==========================================+=================================+
     
	 Also can be used to filter (-f) or extract (-f -R) a list of ids from a file. 


    Usage:
      bio_fileformat_tool.pl [-h] -i <input_file> -a <input_format> -o <output_file> -z <output_format>     
      
    Flags:
      -i   input_file (mandatory)
      -a   input_format (mandatory)
      -o   output_file (mandatory)
      -z   output_format (mandatory)
      -f   a list of id to exclude in the new file (it will ignore everything after the first column)
      -R   use a reverse option with the filter and select only the ids of the filter file
      -h   this help

EOF
exit (1);
}

