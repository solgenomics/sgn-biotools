
use strict;
use Bio::SeqIO;

my $coords_file = shift;
my $fasta_file = shift;
my $size = shift || 3000;

open(my $F, "<", $coords_file) || die "Can't open coords file $coords_file\n";

my %coords;

while (<$F>) {
    chomp;
    my ($original_name, $name, $match_quality, $seq_id, $method, $type, $start_coord, $end_coord, $frame, $strand) = split /\t/;
    print STDERR "Building hash for $original_name, $name\n";
    $coords{$name}->{seq_id} = $seq_id;
    $coords{$name}->{start_coord} = $start_coord;
    $coords{$name}->{end_coord} = $end_coord;
    $coords{$name}->{strand} = $strand;
    $coords{$name}->{original_name} = $original_name;
}


close($F);

my $bio_io = Bio::SeqIO->new( -format => "fasta", -file => $fasta_file );
my $bio_out = Bio::SeqIO->new ( -format => "fasta", -file => ">$fasta_file".".out");

while (my $seq = $bio_io->next_seq()) {

    foreach my $name (keys %coords) {

	if ($coords{$name}->{seq_id} eq $seq->id()) {
	    print STDERR "Processing $name\n";
	    $coords{$name}->{seen} = 1;
	    my $seq_id = $seq->id();
	    # extract sequence
	    my $final_subseq;
	    
	    if ($coords{$name}->{strand} eq "-") {
		my $start =  $coords{$name}->{end_coord}+1;
		my $end = $coords{$name}->{end_coord}+$size;
		my $subseq_string = $seq->subseq($start, $end);
		my $subseq = Bio::Seq->new();
		$subseq->seq($subseq_string);
		$final_subseq = $subseq->revcom();
		$final_subseq->id($name."_".$size."_upstream");
		$final_subseq->desc("original_name: $coords{$name}->{original_name} $coords{$name}->{strand} $coords{$name}->{start_coord} $coords{$name}->{end_coord} excised: $start $end");
	    }
	    else {
		my $start = $coords{$name}->{start_coord} - $size;
		my $end = $coords{$name}->{start_coord} -1;
		my  $subseq_string = $seq->subseq($start, $end);
		my $subseq = Bio::Seq->new();
		$subseq->seq($subseq_string);
		$subseq->id($name."_".$size."_upstream");
		$subseq->desc("original_name: $coords{$name}->{original_name} $coords{$name}->{strand} $coords{$name}->{start_coord} $coords{$name}->{end_coord} excised: $start $end");
		$final_subseq = $subseq;
	    }
	    $bio_out->write_seq($final_subseq);
	    
	}
    }
    

}

foreach my $name (keys %coords) {
    if (! defined($coords{$name}->{seen})) {
	print STDERR "NOT SEEN: $name\n";
    }

}


$bio_io->close();
$bio_out->close();


    
