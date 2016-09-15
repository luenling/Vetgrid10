use strict;
# get variables
my $usage = "Usage:  perl scriptname.pl <infile >";
my $fq = shift or die $usage;
($fq =~ /^-+h/) && die $usage;
# get list of files
my @solexa = '';
my @sanger = '';
my @not_def = '';
open FQ, "<", $fq or die $!;
# read first line
my $line = <FQ>;
my @fields=split("\t",$line);
my $n_pops=(scalar(@fields)-3)/3;
print "Number of populations: ",$n_pops,"\n";
my %pops=();
my $count=0;
while ($line=<FQ>)
  {
  @fields=split("\t",$line);
  for (my $j = 0; $j <=  ($n_pops-1); $j ++){
  check_qual($fields[$j*3+5],$j,\%pops); # divide in chars
  } 
  $count+=1;
  if (scalar(keys(%pops)) >= $n_pops || $count > 10000000) {last;}
}

foreach my $i (sort {$a <=> $b} keys(%pops) ){
print $i,":",$pops{$i},"\n",; 
}


sub check_qual{ # takes a line, the population number and a pointer to the result hash, changes hash 
  my $field=shift;
  if (length($field) <= 1) {return;};
  my $pop=shift; 
  my $hash=shift;
  if (exists(${$hash}{$pop})) {return;};
  my @entries=split(//,$field); # divide in chars
  foreach my $i (@entries){
    my $number = ord($i); # get the number represented by the ascii char
    # check if it is sanger or illumina/solexa, based on the ASCII image at http://en.wikipedia.org/wiki/FASTQ_format#Encoding
    if($number > 74){ # if solexa/illumina and not illumina 1.8
      # "This file is solexa/illumina format\n"; # print result to terminal and die
      ${$hash}{$pop}='sol/il'; last;
    }elsif($number < 59){ # if sanger
      #die "This file is sanger format\n"; # print result to terminal and die
      ${$hash}{$pop}= 'sanger';last;
    }
  }
  return
}
