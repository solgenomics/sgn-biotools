
=head1 NAME

split_fasta_Ns.pl - splits every entry in a fasta file at the Ns and creates a new file with sequence entries representing the split sequences.

=head1 SYNOPSIS

perl split_fasta_Ns.pl fastafile.fa

Creates a new fasta file called fastafile.fa.split.

=head1 AUTHOR

Lukas Mueller (lam87@cornell.edu)

=cut


use strict;

use Bio::SeqIO;

my $file = shift;

my $in = Bio::SeqIO->new( -file => $file, -format => "fasta" );
my $out = Bio::SeqIO->new( -file => ">".$file.".split", -format => "fasta");
while (my $seq = $in->next_seq()) { 
    my $id = $seq->id();
    $id=~s/\|//g;

    my $s = $seq->seq();
    $s=~s/\n//g;
    $s=~s/\r//g;
    my @sub_s = split /N+/, $s;

    my $n = 1;
    foreach my $sub (@sub_s) { 
	print STDERR "Creating new sequence $id.$n...\n";
	my $new_s = Bio::Seq->new();
	$new_s->id("$id.$n");
	$new_s->seq($sub);
        $out->write_seq($new_s);
	$n++;
    }
}
$out->close();
$in->close();

print STDERR "Done.\n";
