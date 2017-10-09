
use strict;

use Bio::SeqIO;

my $file = shift;

my $in = Bio::SeqIO->new( -file => $file, -format => "fasta" );

while (my $seq = $in->next_seq()) { 
    
    my $id = $seq->id();
    $id=~s/\|//g;
    my $out = Bio::SeqIO->new( -file => ">".$file.".$id", -format => "fasta");

    $out->write_seq($seq);

    $out->close();
}

$in->close();
