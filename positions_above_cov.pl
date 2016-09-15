use Getopt::Long;
# use as perl positions_above_cov.pl --mcov "10,100,120,130" < infile.sync > outfile.sync

$mcov=15; # default min coverage
GetOptions ("mcov=s" => \$mcov)   # numeric: minimal coverage of each popualiton required to keep a position
  or die("Error in command line arguments\n");

@covs=split(",",$mcov);
#print @covs;
while (<>){
  $line=$_;
  
  chomp($line);
  @fields=split("\t",$_);
  $coverage=0;
  foreach $j (3..($#fields-$pv)) {
    $sum = 0;
    @entries=split(":",$fields[$j]);
    foreach $i (0 .. 3) {
      $sum +=  $entries[$i];
    }
    ($sum > @covs[$j-3] * 1.1 ) && ($coverage = 1);
  }
  #print $coverage,"\n";
  ($coverage) && (print $_); 
}
