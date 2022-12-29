
use strict;
use Data::Dumper;
use Getopt::Std;
use Bio::SeqIO;

our($opt_f, $opt_e);

my %polys;

getopts('f:e:');

if (! $opt_e) {
    die "-e required. (Tab delimited file with 3 columns: 1 header line, then: polymorphism_id, reference nucleotide, alternate nucleotide)";
}
else {
    open(my $F, "<", $opt_e) || die "Can't open file $opt_e";
    my $header = <$F>;
    
    while (<$F>) {

	chomp;
	
	my ($pm, $ref, $alt) = split /\t/;

	my $id;
	my $coord;
	
	if ($pm !~ m/\_/) {
	    $id = "NONE";
	    $coord = $pm;
	    print STDERR "No sequence identifier ($coord)\n";
	}
	else {
	    
	    ($id, $coord) = split/\_/, $pm;

	    print STDERR "complete id : $id, $coord\n";
	}
	
	$polys{$id}->{$coord}->{ref} = $ref;
	$polys{$id}->{$coord}->{alt} = $alt;
	$polys{$id}->{$coord}->{coord} = $coord;
    }
}

print STDERR "Polymorhpisms: ".Dumper(\%polys);

if (! $opt_f) {
    die "-f required. (Fasta file)";
}

my $seq_count = 0;
my $io = Bio::SeqIO->new( -format=> 'fasta', -file=> $opt_f);

while(my $seq = $io->next_seq()) {
    my $seq_id = $seq->id;
    $seq_count++;
    
    foreach my $id (keys %polys) {

	print STDERR "Processing $id...\n";

	if ($id eq 'NONE' || $id eq $seq_id) {

	    foreach my $c (sort { $a <=> $b } (keys(%{$polys{$id}}))) { 

		print STDERR "COORDINATE: $c\n";
		my $s = $seq->seq();
		
		if ($c+4 > $seq->length()) {
		    $c = $c + ($seq->length() - $c);
		}
		
		my $ref_nucleotide = $seq->subseq($c, $c);
		
		if (uc($ref_nucleotide) ne uc($polys{$id}->{$c}->{ref})) {
		    print STDERR "WARNING! Ref nucleotides do not match ($ref_nucleotide vs. $polys{$id}->{$c}->{ref})\n";
		}
		
		my $left_subseq = $seq->subseq($c-4, $c-1);
		my $right_subseq = $seq->subseq($c+1, $c+4);
		
		print $c.": ".uc($left_subseq). '['.$polys{$id}->{$c}->{ref}."/".$polys{$id}->{$c}->{alt}.']' . uc($right_subseq)."\n";
	    }
	}
    }
    
}

print STDERR "Processed $seq_count sequences. Done.\n";
