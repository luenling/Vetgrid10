#use PDL::GSL::CDF;
#use PDL::Stats::Basic;
use Data::Dumper;
use Math::CDF;
use strict;

sub Welch_Test{
  # takes two population array references (number, sum, and sum of squares of measurements)
  # you can call it either with anonymous arrays: Welch_Test([6,18,58],[8,12,20]) or with two array references  Welch_Test(\@pop1,\@pop2)
  # be aware, however, it will alter the arrays referenced; if you do not want the given arrays to be altered uncomment the lines after my ($pop1,$pop2)=@_;
  # calculates the Welch t test statistic, degrees of freedom, and one sided p Value for the Null-Hypotheses that pop2 has a higher average value than pop1
  # needs the package  Math::CDF from Math-CDF, you can install it with something like this: sudo cpan /Math-CDF/
  # 
  # returns the t value, the degrees of freedom and the P value
  my ($pop1,$pop2)=@_;
  # if you do not want to alter the referenced arrays uncomment the following two lines:
  #$pop1 = [ map{ $_ } @{$pop1} ]; 
  #$pop2 = [ map{ $_ } @{$pop2} ]; 
  #print Dumper($pop1), Dumper($pop2);
  # define array positions for number, mean and variance
  my ($n,$u,$v)=(0,1,2);
  if ($pop1->[$n] <= 1 || $pop2->[$n] <= 1 || ($pop1->[$n] + $pop2->[$n]) < 3 ){
    # not enough reads to say much
    return [0,0,1.0];
  }
  $pop1->[$u]=$pop1->[$u]/$pop1->[$n]; # get $population1 mean
  $pop2->[$u]=$pop2->[$u]/$pop2->[$n]; # get $population2 mean
  #print $pop1->[$u]," ",$pop2->[$u],"\n";
  if ($pop1->[$u] <= $pop2->[$u]){
    # alternative allele has greater tail distance mean than reference
    return [0,0,1.0];
  }
  $pop1->[$v]=($pop1->[$v]-$pop1->[$n]*$pop1->[$u]**2)/($pop1->[$n]-1); # get $pop1 variance
  $pop2->[$v]=($pop2->[$v]-$pop2->[$n]*$pop2->[$u]**2)/($pop2->[$n]-1); # get $pop2 variance
  # calculate t value
  my $t = ($pop1->[$u] - $pop2->[$u])/ sqrt($pop1->[$v]/$pop1->[$n] + $pop2->[$v]/$pop2->[$n]);
  # calculate v (degree of freedoms)
  my $df = ($pop1->[$v]/$pop1->[$n] + $pop2->[$v]/$pop2->[$n])**2/($pop1->[$v]**2/($pop1->[$n]**2*($pop1->[$n]-1)) + $pop2->[$v]**2/($pop2->[$n]**2*($pop2->[$n]-1)));
  # print $t, " ",$df,"\n";
  my $pV = 1.0 - Math::CDF::pt($t,$df),"\n";
  # my $pV =  PDL::GSL::CDF::gsl_cdf_tdist_P($t,$df); 
  #print $pV;
  # my $pV = 1 - gsl_cdf_tdist_P($t,$df); # onesided, Null Hyp: pop2 has higher or equal tail distance
  return [$t,$df,$pV]  
}

# comparison to R:
# a=c(2,3,4,2,3,4)
# b=c(1,2,1,2,1,2,1,2);
# length(a)=6;sum(a)=18; sum(a^2)=58;
# length(b)=8;sum(b)=12; sum(b^)=20
# t = 3.6483, df = 7.645, p-value = 0.003519
my @pop1=(6,18,58);
my @pop2=(8,12,20);
#my $a = Welch_Test([6,18,58],[8,12,20]);
my $a = Welch_Test(\@pop1,\@pop2);
print Dumper($a),Dumper(\@pop1),Dumper(\@pop2);
print "pV: $a->[2]; t: $a->[0]; degrees of freedom: $a->[2]\n";

