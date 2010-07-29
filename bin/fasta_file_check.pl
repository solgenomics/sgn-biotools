#!/usr/bin/perl

=head1 NAME

fasta_file_check.pl - a utility for checking the integrity of fasta formatted files

=head1 SYNOPSIS

perl fasta_file_check.pl [-simvhD] dnafastafile

perl fasta_file_check.pl -p[simvhD] proteinfastafile

Checks the structure of a DNA fasta file (or protein fasta files using the -p option) and prints out various information. For a fast check of a file, use the -s option, which prints a summary report only.

The program checks

=over 4

=item (1)

if letters other than A,T,G,C and N occur in lines that do not start with > (or letters corresponding to the one letter amino acid code +X +* when the protein option (-p) is used). 
Fewer than 80% A,T,G,C and N (or non-amino acid letters) are considered a severe problem 
(such errors will break most parsers etc.).

=item (2)

that there are no empty lines in the files. 

=item (3)

whether duplicate identifiers are present in the file.

=item (4)

whether identifiers exceed a length of 50 characters.

=item (5)

whether there are lines that are over 100000 characters long (the total sequence length of an entry
can be infinitely long - that is not checked. Line length here means a stream of characters until a 
carriage return is encountered. The original fasta definition calls for a line length of 60 characters, but this is not checked as most fasta files don't adhere to this.

=back

Lines with more than 20% non-DNA (or non-protein) characters are counted as severe errors, minor error with less.

=head1 OPTIONS

=over 4

=item -p

Use for files containing protein sequences.

=item -s

Print summary only.

=item -i

Print identifiers only, and summary.

=item -m

Ignore minor problems.

=item -r

Supply a regular expression to use for validating protein/nucleotide sequence.
This can be any regular expression that matches a single character, and
ought to be a character class (e.g., [ACTGXactgx]).  The expression will
be matched case-sensitively, so the user is responsible for ensuring
that both cases are represented, if that is the desired effect.  The default
expressions are

[ACTGactgNn] for nucleotides,
[ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy\*] for protein.

=item -n

print statistics of overall nucleotide frequencies in the file. I would think that this slows the execution considerably, and may increase the memory footprint, if large line lengths are used.

=item -v

Print version of the program.

=item -h

Print help message

=item -D

Display sequence length distribution

=back

=head1 VERSION

 Version 2.5 - Dec 28, 2009

=head2 Version History

 2.0 2002/02/08 Added command line options
 2.1 2002/02/16 Added protein option
                duplicate id checking
                line length checking
                help information
 2.2 2002/09/17 Fixed a bug that reported too many empty lines
 2.3 2004/04/02 Added option D to display sequence length 
                distribution
                Removed option S. Min and max lengths are now 
                always displayed.
 2.4 2007/05/25 Added option to specify a regular expression to 
                match sequence elements against (for unusual 
                program outputs).
 2.5 2009/12/28 Added -n option.

=head1 AUTHOR

 Lukas Mueller

 The Arabidopsis Information Resource. Copyright (c) 2001, 2002. 
 The Carnegie Institution of Washington, Department of Plant 
     Biology. Copyright (c) 2001, 2002.
 Cornell University. Copyright (c) 2004.
 The Boyce Thompson Institute. Copyright (c) 2008-2009.
 The Solanaceae Genomics Network. Copyright (c) 2004-2009.
 All Rights Reversed.

Please report bugs to lam87@cornell.edu.

=cut
    
use strict;
use Pod::Usage;
use Getopt::Std;
use vars qw($opt_s $opt_m $opt_i $opt_v $opt_h $opt_p $opt_r $opt_D $opt_L $opt_n);

# opt_p : protein file.
# opt_s : summary only. Don't generate messages during run.
# opt_n : get stats on nucleotide frequencies
# opt_m : do not report minor problems
# opt_i : identifiers only. Print affected identifiers.
# opt_v : print version info
# opt_h : print help
# opt_D : diplay sequence length distribution
# opt_r : regular expression to use instead of hard-coded nuclotide/amino acid letters.  probably should be a character class.

# constants
#
my $maxnamelength = 50;
my $maxlinelength = 100000; # this is not the length of an entry, just the length of a line 
                            # (of which an entry can have unlimited)
# variable declarations
#
my $lineno = 0;
my $dirtylines = 0;
my $length = 0;  # so first entry doesn't appear to be empty
my $linetoolongcount = 0;
my $sequences_not_checked=0;
my $entryname = "";
my $oldentryname = "";
my $nametoolongcount = 0;
my $zerolengthcount = 0;    # looks at length of a sequence entry
my $zerolinelengthcount =0; # looks at the length of a individual line
my $severeproblemcount = 0;
my $entryProblemCount;
my %entries;
my $severity;
my $duplicatedIds = 0;
my $what;
my $len;
my $good;
my @seqstats;
my $totalentries=0;
my $totallength=0;
my $longest_length=0;
my $shortest_length=0;
my $longest_seq="";
my $shortest_seq="";
my $length_cutoff=0; # used for option L with the length cutoff 
my $longer_seq_count=0; #sequences longer than length_cutoff
my $shorter_seq_count=0; #sequences shorter than length_cutoff
my $last = 0; # some fuctions of process_entry should not be preformed when it is called a final time - used as marker

my $n_ref = {}; # hashref with counts of different symbols

getopts('psmnivhr:DL:');

if ($opt_v) { print_version(); exit(); }

if ($opt_h) {  pod2usage(-verbose => 1); }
 

my $filename = shift;

open (F, "<$filename") || die "Can't open $filename. Execution terminated. Please try again.\n";

print STDERR "Checking file... \n";
							      
my $line;
							      
while (<F>) {
   chomp;
   $line = $_;
   check_line($line, $n_ref);
   
}

$last = 1; ## marks this sub call as the final call
process_entry($line);

close (F);

# quit with option L without giving the summary, as this would mess up the tab delimited output 
# when it is used in sorting etc.
if ($opt_L) { 
    print STDERR "# >= $opt_L: $longer_seq_count\n";
    print STDERR "# <  $opt_L: $shorter_seq_count\n";
    exit(0); 
}

print "\nSummary\n*******\n\n";
print "This file contained the following errors:\n";
if ($opt_m) {
    print "- [Minor errors ignored (option -m)]\n";
}
else { 
    print "- $dirtylines lines with minor sequence problems\n";
}
print "- $severeproblemcount lines with major sequence problems\n";
print "- $nametoolongcount entries where the seq name is too long\n";
print "- $zerolengthcount entries that lack sequence info.\n";
print "- $linetoolongcount lines that are over $maxlinelength long.\n";
print "- $zerolinelengthcount lines that are emtpy [may be trailing]\n";
print "- $duplicatedIds duplicated Ids.\n";
print "- $sequences_not_checked sequences were not checked [too long]\n\n";
print "- Total length: $totallength\n";
print "- Sequence entries: $totalentries\n";
print "- Average Sequence length: "; printf "%10.2f", ($totallength/$totalentries);
print "\n";
print "- Longest Sequence: $longest_length\n";
print "- Shortest Sequence: $shortest_length\n";
print "\n";

if ($opt_n) { 
    my $total_nucs = $n_ref->{C} + $n_ref->{G} + $n_ref->{A} + $n_ref->{T} + $n_ref->{N};

    print "Residue Stats:\n";
    print "==============\n";

    foreach my $k (keys %$n_ref) { 
	print $k . " " . $n_ref->{$k}." ".(sprintf "%5.3f", ($n_ref->{$k}/$total_nucs)*100)."\%\n";
    }


    if (! $opt_p) { 

	print "Overall GC content: ". (sprintf "%5.3f", ((($n_ref->{C}+$n_ref->{G})/$total_nucs * 100))). " \n";
    }
}

my $exit_code = 0;

if ($severeproblemcount || $zerolengthcount) { print "This file should be checked and repaired before further use. \n\n"; $exit_code=-1;  }
elsif ((!$severeproblemcount) && !($nametoolongcount) && (!$zerolengthcount)) { print "This file should be usable.\n\n"; }

if ($opt_D) {
  print "Sequence Size Distribution:\n\n";
  print "  Size   Count\n";
  my $i=0;
  my @compressed =();
  my %label = ();
  my $divs = 20;
  my $increment = @seqstats/$divs;

  for (my $i=0; $i<@seqstats; $i++) {
      $compressed[$i/$increment]+=$seqstats[$i];
  }
  my $max = 0;
  foreach my $i (@compressed) {
      if ($max < $i) { $max= $i; }
  }
  my $maxwidth = 60;
  my $factor=$maxwidth/$max;

#  print "max=$max. factor=$factor\n";
  for (my $i=0; $i<@compressed; $i++) {
      my $size=$i*$increment;
      printf "%6d",  ($size); print " ["; printf "%6d", "$compressed[$i]"; print "]";
      print " ";
      for (my $n=0; $n<($compressed[$i] * $factor); $n++) {
	  print "*";
      }
      print "\n";
  }
}
 

 
      
print "\nDone.\n\n";

exit($exit_code);

sub checkseq {     
    $what = $_[0];
    $len = length($what);
    if ($len > 32000) { return 3; }
    my $re = $opt_r ? $opt_r : "[ACTGactgNnXx]";
    if ($what =~ /${re}{$len}/) {
	return 0;
    }
    else {
	$_ =  $what;
	$good = tr/$re/x/;
	if ($good/$len > 0.8) { return 1;  }
	else { return 2; }
    }
}


sub checkproteinseq {
    $what = $_[0];
    $len = length($what);
    if ($len > 32000) { return 3; }
    my $re = $opt_r ? $opt_r : "[ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy\*]";
    if ($what =~ /${re}{$len}/)
    {
	return 0;
    }
    else
    {
	$_ =  $what;
	$good = tr/$re/x/;
	if ($good/$len > 0.8) { return 1;  }
	else { return 2; }
    }
}

sub message {
    my $message = shift;
    # if 'summary only' is not on, we print everything
    if (!$opt_s) {
	# check if 'identifiers only' (opt_i) option is set
	if ($opt_i) {
	    # print the identifier only once if it contains several problems
	    if ($entryProblemCount <= 1)  { 
		print "$entryname ($severity)\n";
	    }
	}
	else { print $message; }
    }   
}
    
sub print_version {
    print "\nfasta_file_check Version 2.4 May 25, 2007\n";
    print "Copyright (c) 2001-2004 by Lukas Mueller,\n";
    print "The Arabidopsis Information Resource - TAIR, and \n";
    print "The Carnegie Institution of Washington, Department of Plant Biology.\n\n";
}

sub print_help {
    print "\nfasta_file_check\n\nHelp\n====\n\nOptions:\n";
    print "-p : protein file.\n";
    print "-s : summary only. Don't generate messages during run.\n";
    print "-m : do not report minor problems.\n";
    print "-i : identifiers only. Print sequence identifiers of error containing entries.\n";
    print "-r : supply a regular expression to use for validating sequences\n";
    print "-v : print version info.\n";
    print "-h : print help (these lines, actually).\n";
    print "-D : display sequence length distribution.\n";
    print "-L : print out a tab delimited file with the lengths of the sequences and their ids (no other output).\n\n";


}

sub hashValueAscendingNum {
   $a <=> $b;
}

sub check_line {
    my $line = shift;
    my $n_ref = shift;

    $lineno++;
    if (length($line) == 0) {
	message("Problem at line $lineno: Line has zero length!\n");
       $zerolinelengthcount++;
       $severeproblemcount++;
    }
    if (($line =~ /^>/)) { 
       process_entry($line);
   }
   else { 
       $length += length($line);   
       if ($opt_p) { $severity = checkproteinseq($line); }
       else {
	   $severity = checkseq($line);
       }
       if ($severity >0)  { 
	   if (($severity ==1)) {
	       # do we show minor problems? check m option
	       if (!$opt_m) {
		   # if minor problem and option 'ignore minor problem' not set, print it 
		   message("Problem at line \# $lineno :\n$line\n\n");
		   $dirtylines++; 
		   $entryProblemCount++;
	       }
	   }
	   elsif ($severity ==2) { 
	       message("Severe problem at line \#$lineno: \n$line\n\n");
	       $severeproblemcount ++;
	       $entryProblemCount++;
	   }
	   elsif ($severity == 3) {
	       message("Sequence entry too long to be checked: $lineno $entryname\n");
	       $sequences_not_checked++;
	   }
	   else { 
	       message ("Unknown problem at line \# $lineno. Severity: $severity.\n\n"); 
	   }
       }
       else {
	   ## do nothing because sequence is fine
       } 
       $severity = 0;
       if ($opt_n) { 
	   get_nuc_stats($line, $n_ref);
       }
       
   }       
}

sub process_entry {

    my $line = shift;
    if ($last == 0) { ## won't add to this total when sub is called an extra time
	$totalentries++;
    }
     $totallength+=$length; ## this IS needed to add the last sequence length to total
     if ($length>$longest_length) { $longest_length = $length; }
     if (($length<$shortest_length) || ($shortest_length==0)) { $shortest_length = $length; }
     if (($length ==0) && ($lineno!=1)) { 
	 $entryProblemCount++;
	 message ("Problem at entry $entryname (line# $lineno) No sequence found!!! \n\n"); 
	 $zerolengthcount++; 
	 $severeproblemcount++;
     }
     if (length($line)>$maxlinelength) { 
	 message("Problem at line $lineno: Line too long (".length($line).">$maxlinelength)\n"); 
	 $linetoolongcount++;
     }
     if ($opt_D) {
	 if ($length != 0) { $seqstats[$length]++; } ## won't add initialized $length value to array
     }
     if ($opt_L) {
	 print "$entryname\t$length\n";
	 if ($length >= $opt_L) { $longer_seq_count++; }
	 if ($length < $opt_L) { $shorter_seq_count++; }
     }
     
     $length = 0; 
     $entryProblemCount = 0;
    if ($last == 0 ) {
	$oldentryname = $entryname;
	$entryname = $line;
	$entryname =~ s/\>(.*?)\s(.*)/$1/;

	# check if entries are duplicated.
	if (exists ($entries{$entryname})) { 
	    message("$entryname is duplicated\n");
	    $duplicatedIds++;
	}
	else {
	    $entries{$entryname}++; 
	}
	if ( length($entryname) > $maxnamelength) { 
	    message("Problem at entry $entryname (line# $lineno): Name too long!\n"); $nametoolongcount++; 
	}
    }

}

sub get_nuc_stats { 
    my $line = shift;
    my $n_ref = shift;

    my @line = split //, $line;
    foreach my $c (@line) { 
	$n_ref->{uc($c)}++;
    }
}




