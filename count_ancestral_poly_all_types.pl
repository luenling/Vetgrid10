use Getopt::Long;
$winsize=5e5;
GetOptions ("w=i" => \$winsize)   # numeric
  or die("Error in command line arguments\n");
$win=1;
$chr=0;
while (<>){
  chomp();
  @fields=split("\t");
  @mono=(0,0,0); # monomorph sim 0: not conserved, 1: major conserved, 2: minor conserved 
  @poly=(0,0,0,0,0); # polymorph sim: 0: ncs, 1: maj cons 2: min cons 3:  4: 
  $coverage = 1; 
  foreach $j (3..$#fields) {
    $sum = 0;
    @entries=split(":",$fields[$j]);
    foreach $i (0 .. 3) {
      $sums[$i] += $entries[$i];
      $sum +=  $entries[$i];
    }
    ($sum < $mcov) && ($coverage = 0);
   }
  @sort_sums = sort {$b <=> $a} @sums;
  @sort_sums = sort {$b <=> $a} @sums;
  if ($sort_sums[1] >= 5  && $coverage) {print $_}; 
}
