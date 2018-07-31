
=head1 NAME

stitch.pl - stitch together fasta files from ordered pieces

=head1 DESCRIPTION

perl stitch.pl [fasta file] [match file blast m8]

Output will be in files named after chromosomes: stitched[chr-nr].fa


=head1 AUTHOR

Lukas Mueller <lam87@cornell.edu>

=cut


use strict;
use Bio::SeqIO;
use File::Slurp qw | read_file |;

my $seq_file = shift;
my $match_file = shift;

my $in = Bio::SeqIO->new( -format => "fasta", -file => $seq_file);

my %seqs;

my $count = 0;
while (my $seq = $in->next_seq()) { 
    print STDERR "\rReading sequence file... $count";
    $seqs{$seq->id()} = $seq->seq();
    $count++;
}
print STDERR " Done.\n";

my $n_seq = 'N' x 100;

print STDERR "Sequence spacer is: $n_seq\n";

my $previous_query = "";

my @matches;

if (-e $match_file.".parsed") { 
    print STDERR "Reading cached blast data... ";
    my @lines = read_file($match_file.".parsed");
    chomp(@lines);
    foreach my $line (@lines) { 
	push @matches, [ split /\t/, $line ];
    }
}

else { 
    open(my $F, "<", $match_file) || die "Can't open match file $match_file\n";
    
    
    
    $count = 0;
    while (<$F>) { 
	chomp;
	my ($query, $subject, $percent, $aln_length, $a, $b, $q_start, $q_end, $s_start, $s_end, $evalue, $score) = split /\t/;
	
	if ($previous_query ne $query) { 
	    my $orientation = '+';
	    if ($s_end < $s_start) { $orientation = '-'; }
	    
	    my $s_start_coord = 0;
	    if ($orientation eq '+') { 
		$s_start_coord = $s_start - $q_start;
	    }
	    else { 
		$s_start_coord = $s_start - (length($seqs{$query}) - $q_end);
	    }
	    
	    $subject=~s/lcl\|Contig(\d+).*/$1/g;
	    push @matches, [ $query, $subject, $s_start_coord, $orientation, $q_start, $q_end, $s_start, $s_end, $percent, $aln_length, $evalue, $score ];
	    $previous_query = $query;
	}
	if ($count % 1000000 == 0) { print STDERR "\rParsing blast report... ".int($count/1000000)."M lines ";}
	$count++;
    }
}

print STDERR "Done.\n";

print STDERR "Sorting matches... ";
my @sorted_matches = sort { if ($a->[1] == $b->[1]) { $a->[2] <=> $b->[2] } else {  $a->[1] <=> $b->[1]  } } @matches;
print STDERR "Done.\n";

print STDERR "Write sorted assembly file... ";

if (!-e $match_file.".parsed") { 
    my $G;
    open($G, ">", $match_file.".parsed") || die "Can't write file $match_file.parsed";
    
    foreach my $line (@sorted_matches) { 
	print $G join( "\t", @$line);
	print $G "\n";
    }
}

my $previous_subject = "";
my $H;
my $open_file_handle =0;
foreach my $line (@sorted_matches) { 
    my ($query, $subject, $start, $orientation) = @$line;

    print STDERR "$subject, $previous_subject, $open_file_handle\n";    
    if ($subject ne $previous_subject) { 

	if ($open_file_handle) { close($H); }
	open($H, ">", "stitched_chr$subject.fa") || die "Can't open chr file\n";
	$open_file_handle=1;
	print $H ">seq$subject\n";
	print "\n>seq$subject\n";
    }

    my $seq = $seqs{$query};
    
    if ($orientation eq "-") { 
	print STDERR "Using reverse complement seq for $query...\n";
	my $so = Bio::Seq->new();
	$so->seq($seqs{$query});
	my $revcom = $so->revcom();
	$seq = $revcom->seq();
    }
    else { 
	print STDERR "Forward seq for $query...\n";
    }

    if (! $open_file_handle) { die "ouch!\n"; }

    print $seq;
    print $n_seq;
    print $H $seq;
    print $H $n_seq;
    
    $previous_subject =  $subject;
}
close($H);

print STDERR "Done.\n";

    
