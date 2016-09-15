#!/usr/bin/perl
# author: Lukas Endler Tue Dec 18 14:00:07 CET 2012
# takes an inputstream in syncfile format and calculates average coverage and stdv chrosome wise and overall for each population 
use warnings;
use strict;
#use Data::Dumper;
#use Getopt::Long;
#use Pod::Usage;


my $threshold=5;
my @chrom_levs = qw( Sum SSum Zeros Length Mean Stdv );
my %chromnames=();
my @chroms=();
my $chromcount=0;
my $count = 0;
my $pops = 0;
my $spaces = qr/\s+/o;
my $lastchrom = "";
print STDERR "reading sync file:\n"; 
while (<>){
  $count += 1;
  #my $line = $_;
  chomp($_);
  #my @fields = split(/\s+/, $line);
  my @fields = split($spaces, $_);
  my $c_name = $fields[0];
  if ($c_name ne $lastchrom) {
    $lastchrom = $c_name;
    print STDERR "Chrom $c_name : 1 \n";
    unless (exists $chromnames{$c_name}){
      $chromcount = push(@chroms,[ $c_name, [ ] ]);
      $chromnames{$c_name} = $chromcount - 1;
      $chroms[$chromcount -1][1] = [ undef,undef,undef,0,0,0 ];
      # Sum=0 SSum=1 Zeros=2 Length=3 Mean=4 Stdv=5
    }
  }

  my $cindx = $chromnames{$c_name};
  $chroms[$cindx][1][3] += 1;
  unless ($chroms[$cindx][1][3] % 1e5) {print STDERR "Chrom $c_name : ",$chroms[$cindx][1][3],"\n"}; 
  $pops = $#fields -3;
  unless (defined $chroms[$cindx][1][0]) {  
    for (0..2) { $chroms[$cindx][1][$_] =  [ map {0} (0..$pops) ] }; 
  }  
  for my $i (3..$#fields) {
    my $val = 0;
    $val += $_ for (split(":",$fields[$i]));
   # for my $j (split(":",$fields[$i])) { $val += $j};
    $chroms[$cindx][1][0][$i-3] += $val;
    $chroms[$cindx][1][1][$i-3] += $val*$val;
  #  if ($val == 0) {$chroms{$c_name}{'Zeros'}->[$i-3] += 1 };
  }
 # ($count == 100) && last;
} 

print STDERR "finished reading, calculating coverage\n"; 
for my $cindx (0..$#chroms){
  $chroms[$cindx][1][4] = mean($chroms[$cindx][1][0],$chroms[$cindx][1][3]);
#  print STDERR "$pops length $chroms[$cindx][1][3] mean  $chroms[$cindx][0] : @{$chroms[$cindx][1][4]}  \n";
#  print STDERR Dumper($chroms[$cindx][1][0]);
  $chroms[$cindx][1][5] = stdv($chroms[$cindx][1][0],$chroms[$cindx][1][1],$chroms[$cindx][1][3]);
}


print "chrom\t",join("\t",map(sprintf("Pop%i",$_+1), (0..$pops))),"\tLength\n"; 

for my $cindx (0..$#chroms){
  print $chroms[$cindx][0]; 
  for my $j (0..$#{$chroms[$cindx][1][4]}){
	printf("\t%.1d %.1d",$chroms[$cindx][1][4][$j],$chroms[$cindx][1][5][$j])    
    }   
  printf("\t%i\n",$chroms[$cindx][1][3]);
}
## Calculate totals
## tot. length
my $overall_length = 0;
$overall_length += $chroms[$_][1][3] for ( 0..$#chroms );
## total mean
my @overall_mean = map {my $i=0; for my $cindx (0..$#chroms) {$i += $chroms[$cindx][1][4][$_]*$chroms[$cindx][1][3]/$overall_length }; $i} (0..$pops);
## total stdv
my @overall_stdv =  map {my $i=0; for my $cindx (0..$#chroms)  {$i += $chroms[$cindx][1][1][$_]}; $i = sqrt($i/$overall_length - $overall_mean[$_]**2)} (0..$pops);
## Print Totals
print "Total" ;
for my $j (0..$pops){
  printf("\t%.1d %.1d",$overall_mean[$j],$overall_stdv[$j])    
}
printf("\t%i\n",$overall_length);

#######################
### end of main
#######################
sub mean {
  (my $sums , my $length) = @_;
  my @means = map { $_ / $length }  @{$sums};
  return \@means ; 
}

sub stdv {
  (my $sums, my $ssums, my $length ) = @_;
  my @stdvs = map {sqrt ( (${$ssums}[$_] - (${$sums}[$_]**2)/$length)/$length  )} (0..$#{$sums});
  return \@stdvs; 
}

