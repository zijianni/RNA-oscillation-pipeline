#!/usr/bin/perl

#the program expects to receive data in a tab-delimited file, with 1 column containing labels for the time steps, and each additional column containing the data for one experimental condition over all the time steps.
#any row in the file that begins with a "#" character will be printed out exactly as it is read in (e.g., a header row). ALL such rows will be printed out before any output data (i.e., if your input file has a "footer" row at the end of the luminescence data, it will NOT be printed out at the end of the output file. instead it will be printed out before the first row of output data.
#./bandpass.pl < my_cool_data.tab 45 3 > my_cool_filtered_data.tab


use strict;
use List::Util qw( min max );

mainSub();
exit;

sub mainSub {
  my (@args)=(@ARGV);
  my $highPassWindow = $args[0] || 45;
  my $lowPassWindow = $args[1] || 3;

  my @lpWeight=();
  my $lpVariance = findVar($lowPassWindow);
  for my $nn (0..$lowPassWindow) {
     $lpWeight[$nn]=normal(0,$lpVariance,$nn);
  }
#print STDERR $lpVariance."\n";
#print STDERR join("\n",@lpWeight,"");
  my $nTSteps=0;
  my @k=();
  my @tStep=();
  while (<STDIN>) {
    if (! m/^#/) {
      chomp;
      my ($tStep,@f)=split /\t/;
      $tStep[$nTSteps]=$tStep;
      my $i = 0;
      for my $f (@f){
        $k[$nTSteps][$i++] = $f;
      }
      $nTSteps++;
    }
    else {print}
  }

  my $nSamples = scalar @{$k[0]};
#print STDERR join("\t", "nSamples", $nSamples, "\n");
  my @kk=();
  my @kkk=();
  for (my $iSample = 0; $iSample < $nSamples; $iSample++) {
    for (my $t = 0; $t < $nTSteps; $t++) {
      my $t0 = max(($t-$highPassWindow), 0);
      my $t1 = min( ($t+$highPassWindow),($nTSteps-1));
      my $kAverage = 0;
      for my $tt ($t0..$t1) {
        $kAverage += $k[$tt][$iSample];
      }
#print STDERR join("\t", $iSample, $t, $k[$t][$iSample], $kAverage, $t1, $t0, $kAverage/($t1-$t0+1),":\n");
      $kAverage/=($t1-$t0+1);
#      $kAverage/=(2 * $highPassWindow)+1;
      $kk[$t][$iSample]=$k[$t][$iSample]-($kAverage);
#print STDERR join("\t", $iSample, $t, $k[$t][$iSample], $kAverage, $kk[$t][$iSample],":\n");
    }
#    for (my $t = 0; $t < 0; $t++) {
    for (my $t = 0; $t < $nTSteps; $t++) {
      my $nMinus= ($t >= $lowPassWindow) ? $lowPassWindow : $t;
      my $nPlus = (($t+$lowPassWindow) < $nTSteps) ? $lowPassWindow : ($nTSteps-$t-1);
      my $kWtdSum = 0;
      my $weightSum  = 0;
      if ($nMinus > 0) {
        for my $tMinus (1..$nMinus) {
          $kWtdSum += $kk[$t-$tMinus][$iSample] * $lpWeight[$tMinus];
          $weightSum += $lpWeight[$tMinus];
        }
      }
      for my $tPlus (0..$nPlus) {
        $kWtdSum += $kk[$t+$tPlus][$iSample] * $lpWeight[$tPlus];
          $weightSum += $lpWeight[$tPlus];
      }
#      $kAverage /= (1+$t0-$t1);
      $kkk[$t][$iSample]=sprintf("%.1f",($kWtdSum/$weightSum));
    }
  }

  my $nExperiments = scalar @{$kk[0]};
  for (my $t = 0; $t < $nTSteps; $t++) {
    print join("\t", $tStep[$t], @{$kkk[$t]})."\n";
  }
  
}

sub findVar {
  my $n=shift @_;

  my $done  = 0;
  my $upper  = $n**2;
  my $var = $upper / 2;
  my $lower = 0.0;
  my $sum = 0.0;
  my $preVar = $upper;
  while (! $done) {
    $sum = 0.0;
    for my $i (-$n..$n) {
      my $nor= normal(0,$var,$i);
      $sum+=$nor;
    }
#print join(":\t:", $done, $n, $lower, $var, $upper,  ($upper-$lower), $sum, normal(0,$var,10),0.01/($n+$n+1),"","\n");
    if ($sum >= 0.999) {
        $lower = $preVar = $var;
        $var = ($var + $upper)/2;
    }
    else {
        $upper = $preVar = $var;
        $var = ($var + $lower)/2;
    }
    if (0.999999999 < ($lower/$upper)){$done=1}
  }
#print join(":\t:", $done, $lower, $var, $upper, $sum, normal(0,$var,10),"","\n");

  return $var;
}

sub normal { 
  my ($mean, $var, $x) = (@_);
  return undef unless $var > 0; 
  return exp(-($x - $mean)**2/(2 * $var))/sqrt(8 * atan2(1,1) * $var);
} 
