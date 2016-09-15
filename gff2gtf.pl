##!/usr/bin/perl
 
use strict;
use warnings;
use Data::Dumper;
 
use File::Basename;
use Bio::Tools::GFF;
 
my $inFile = shift;
my ($name, $path, $suffix) = fileparse($inFile, qr/\.gff/);
my $outFile = $path . $name . ".gtf";
 
my $inGFF = Bio::Tools::GFF->new( '-file' => "$inFile",
 '-format' => 'GFF',
 '-version' => 3 );
my $outGTF = Bio::Tools::GFF->new( '-file' => ">$outFile",
 '-format' => 'GFF',
 '-version' => 2.5);
#my $gffstr = $inGFF->_gff2_string;
while (my $feature = $inGFF->next_feature() ) {
 
print $inGFF->_gff2_string($feature),"\n";
#$outGTF->write_feature($feature);
 
}
