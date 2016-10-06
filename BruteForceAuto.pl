use strict;
use warnings;

use LWP::Simple;
use FileHandle;
use Data::Dumper qw(Dumper);
use Scalar::Util qw(looks_like_number);
use List::Util qw[min max];

######################## Command Line Arguments ################################

# quit unless we have the correct number of command-line args
my $num_args = $#ARGV + 1;
if ($num_args != 5) {
    print "\nUsage: PerlTest.pl RD.out_file_name Epsilon MinPts CutoffStart CutoffEnd\n\n";
    exit;
}

my $inputFile = $ARGV[0];
my $Epsilon = $ARGV[1];
my $MinPts = $ARGV[2];
my $CutoffStart = $ARGV[3];
my $CutoffEnd = $ARGV[4];

##########################  Reading from RD.out  #############################

my @InitialSet;

open(IN, "<$inputFile") || die "Can't open $inputFile: $!";
while (my $line = <IN>) {
	chomp $line;
	my @tabs2 = split(/\t/,$line);
	push @InitialSet, [$tabs2[0],$tabs2[1]];
}

###############################################################################

my @ClusterArray;

for (my $Cutoff = $CutoffStart; $Cutoff <=  $CutoffEnd; $Cutoff += 0.5) {
	my %Clusters;
	my $n = 0;
	my $counter = 0;
	my $s = 0;
	while ($n < scalar @InitialSet) {
		if ($InitialSet[$n][1] < $Cutoff) {
			$counter++;
			if ($counter >= $MinPts) {
				$Clusters{$s} = $n;
			}
		$n++;
		}
		else {
			$s = $n;
			$counter = 0;
			$n++;
		}
	}

	my @OrderedS = sort { $Clusters{$a} <=> $Clusters{$b} } keys %Clusters;
	my @OrderedE = @Clusters{@OrderedS};
	
	foreach my $start (@OrderedS) {
		my $ClusTot = 0;
		for (my $t = $start; $t <= $Clusters{$start}; $t++) {
			$ClusTot = $ClusTot + $InitialSet[$t][1];
		}
		my $ClusAvg = $ClusTot/($Clusters{$start}-$start);
		push @ClusterArray, [$start,$Clusters{$start},$ClusAvg];
	}
}

# print "Cluster Array:\n";
# print Dumper \@ClusterArray;


########################  Writing to file  ##############################

my $OutFile1 = "$inputFile.BruteForceAutoClusters.out";
open (OUT, ">$OutFile1");

for (my $i = 0; $i < scalar @ClusterArray; $i++) {
	print OUT "$ClusterArray[$i][0]\t$ClusterArray[$i][1]\t$ClusterArray[$i][2]\n";
}
close (OUT);

############################  Plotting  #################################

system ("Rscript BruteForceAutoClustersLines.R $inputFile $inputFile.BruteForceAutoClusters.out BruteForceAutoClusters.$inputFile.pdf $Epsilon $MinPts");

print "Done.\n";