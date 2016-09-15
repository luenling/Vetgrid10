#!/usr/bin/perl
#
# VMA2FASTA R2
# 
# Converts nfo format to a fasta file
#
# Copyright (c) 2005 DPGP                                                    
#                                                                              
# This piece of software was created as part of the Drosophila  Population 
# Genomics Project opensource agreement.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell    
# copies of the Software, and to permit persons to whom the Software is        
# furnished to do so, subject to the following conditions:                     
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#

use Getopt::Std;
###############################################################
# Globals
###############################################################
my $gFASTA;      # FASTA w/o indels
my $gFASTAI;     # FASTA with indels
my $gQUAL;
my $gQUALI;


###############################################################
# Main
###############################################################

getopt("iot");


###############################################################
# Verify and process command line parameters
###############################################################
unless ( -e $opt_i && $opt_o ) {
        usage();
}

$opt_t = 15 unless ( defined $opt_t && $opt_t >= 0 && $opt_t <= 99 );

###############################################################
# Open input file
###############################################################
open IN, $opt_i;
open OUT, ">$opt_o.fa" || die;
open OUTI, ">$opt_o.ifa" || die;
open QUAL, ">$opt_o.qual" || die;
open QUALI, ">$opt_o.iqual" || die;

my $i = 0; # last printed genome position
while (<IN>) {
   chomp;
   my $row = $_;
   @cols = split /\t/m, $row;
   ##############
   # Parse Line
   ##############
   my $rp = $cols[0];
   my $bp = $cols[2];
   my $bpq = int $cols[3];
   
   
   ##############
   # Insert Case
   ##############   
   if ( $rp eq 'INSERT' ) {
   	if ( $bpq >= $opt_t ) {
   		$gFASTAI .= lc $bp;
   		$gQUALI .= ( $bpq > 9 ) ? " $bpq" : "  $bpq";
   	}
   	else {
   		$gFASTAI .= 'n';
   		$gQUALI .= '  0';
   	}
   }
   ##############
   # Match Case
   ##############    
   else {
    ##############
    # Catch Up
    ##############
    my $c = int $rp; # Current Genome Position
    while ( $i < $c - 1 ) {
    	$gFASTAI .= 'N';
    	$gFASTA .= 'N';
    	$gQUAL .= '  0';
    	$gQUALI .= '  0';
    	$i++;
   	}
	if ( $bpq >= $opt_t ) {
   		$gFASTAI .= uc $bp;
   		$gQUALI .= ( $bpq > 9 ) ? " $bpq" : "  $bpq";
   		if ( $bp ne '-' ) {
   			$gFASTA .= uc $bp;
   			$gQUAL .= ( $bpq > 9 ) ? " $bpq" : "  $bpq";
   		}
   		else {
   			$gFASTA .= 'N';
   			$gQUAL .= '  0';
   		}
   	}
   	else {
   		$gFASTAI .= 'N';
   		$gQUALI .= '  0';
   		$gFASTA .= 'N';
   		$gQUAL .= '  0';
   	}
   	$i = $c;
   	if ( $i != length $gFASTA ) {
   		die "$i, ", length $gFASTA, " $row\n";
   	}
   }
}


$gFASTA =~ s/(.{60})/$1\n/g;
chomp $gFASTA;
$gQUAL =~ s/(.{180})/$1\n/g;
chomp $gQUAL;
$gFASTAI =~ s/(.{60})/$1\n/g;
chomp $gFASTAI;
$gQUALI =~ s/(.{180})/$1\n/g;
chomp $gQUALI;


print OUT ">$opt_i threshold = $opt_t\n";
print OUT "$gFASTA\n";
print QUAL ">$opt_i threshold = $opt_t\n";
print QUAL "$gQUAL\n";
print OUTI ">$opt_i threshold = $opt_t\n";
print OUTI "$gFASTAI\n";
print QUALI ">$opt_i threshold = $opt_t\n";
print QUALI "$gQUALI\n";

###############################################################
# Sub: Usage
############################################################### 
sub usage {
 my $u = <<END;
VMA2FASTA    -i input vertical multiple alignment .vma
             -o output prefix for .fa and .qual files
             -t score threshold [default 15]
END

 print $u;

 exit(1);
 
}
