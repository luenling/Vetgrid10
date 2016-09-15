use Getopt::Long;
# use as perl polymorph_sites_from_sync.pl --mc 5 --mcov 10 < infile.sync > outfile_mc_5_mcov_10.sync

$mc=5; # default min count
$mcov=15; # default min coverage
$pv=0; # PV in file (cmhout file)
$perc=0; # minimal minor allele freq in each population
GetOptions ("mc=i" => \$mc,    # numeric: minimal overall count of minor allele required to keep a position
	   "mcov=i" => \$mcov,   # numeric: minimal coverage of each popualiton required to keep a position
	   "perc=f" => \$perc,   # float < 1: minimal minor allele freq that needs to be reachein in at least one population required to keep a position
	   "p" => \$pv)   # logic: is pV in sync
  or die("Error in command line arguments\n");

while (<>){
  $line=$_;
  chomp($line);
  @fields=split("\t",$_);
  @sums=(0,0,0,0);
  $coverage = 1;
  $perc_reached=0;
  foreach $j (3..($#fields-$pv)) {
    $sum = 0;
    @entries=split(":",$fields[$j]);
    if ( $perc )  {
	      ($maj,$min,undef) = sort {$b <=> $a} @entries;
	      ((1.0*$min)/($maj+$min)) >= $perc && {$perc_reached = 1}; 
	     }
    foreach $i (0 .. 3) {
      $sums[$i] += $entries[$i];
      $sum +=  $entries[$i];
    }
    ($sum < $mcov) && ($coverage = 0);
   }
  @sort_sums = sort {$b <=> $a} @sums;
  if ($sort_sums[1] >= $mc  && $coverage && ( ! $perc || $perc_reached )) {print $_}; 
}
