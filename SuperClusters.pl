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
if ($num_args != 4) {
    print "\nUsage: SuperClusters.pl RD.out_file_name Epsilon MinPts Eta\n\n";
    exit;
}

my $inputFile = $ARGV[0];
my $Epsilon = $ARGV[1];
my $MinPts = $ARGV[2];
my $Eta = $ARGV[3];
#my $CutoffEnd = $ARGV[4];

##########################  Reading from RD.out  #############################

my @InitialSet;

open(IN, "<$inputFile") || die "Can't open $inputFile: $!";
while (my $line = <IN>) {
	chomp $line;
	my @tabs2 = split(/\t/,$line);
	push @InitialSet, [$tabs2[0],$tabs2[1]];
}

###############################################################################

my %Clusters;
my @ClusterArray;


# Identify super clusters:

my $scs = 0; # super cluster start
for (my $i = 1; $i < scalar @InitialSet; $i++) {
	#print "i=$i\n";
	if ( $InitialSet[$i][1] == 10 ) {
		#print "RD($i)=10\n";

		if ( $InitialSet[$i-1][1] == 10 ) {
			$scs = $i;
		}
		else {
			$Clusters{SuperClusters}{$scs} = $i;
			$scs = $i;
		}		
	}
}
print Dumper \%Clusters;

# Go inside each super cluster and start scanning from smallest r(x):

foreach my $superC ( keys $Clusters{SuperClusters} ) {
	#print "Start of SC = $superC, $InitialSet[$superC][0]\n";
	my $MinRDat = $superC; # Lowest r(x) is at this x value.
	for (my $i = $superC; $i <= $Clusters{SuperClusters}{$superC} ; $i++) {
		if ($InitialSet[$i][1] < $InitialSet[$MinRDat][1]) {
 			$MinRDat = $i;
 		}
	}
	#print "\tMinRD=$InitialSet[$MinRDat][0],$InitialSet[$MinRDat][1]\n";
	my $MaxRD = 10;
	for (my $i = $superC+2; $i < $Clusters{SuperClusters}{$superC} ; $i++) {
		if ($InitialSet[$i-1][1] < $InitialSet[$i+1][1]) {
 			$MaxRD = $InitialSet[$i][1];
 		}
	}

	# Start scanning from epsilon prime= minimum r(x) to maximum r(x):

	my $nsubc = 0; # number of sub-clusters
	my $nsubcPre = 0; # previous number of sub-clusters
	my $TempSubClusterRef;

	for (my $Cutoff = $InitialSet[$MinRDat][1] + 0.1; $Cutoff <=  $MaxRD + 0.5 ; $Cutoff += 0.1) {
		my $n = $superC;
		my $counter = 0;
		my $s = $superC;
		while ($n <= $Clusters{SuperClusters}{$superC}) {
			if ($InitialSet[$n][1] < $Cutoff) {
				$counter++;
				if ($counter >= $MinPts) {
					$Clusters{SubClusters}{$superC}{$s} = $n;
				}
			$n++;
			}
			else {
				$s = $n;
				$counter = 0;
				$n++;
			}
		}
		
		if (exists $Clusters{SubClusters}{$superC}) {  # if at leat one sub-cluster is found in the corresponding super cluster:
			$nsubc = scalar keys $Clusters{SubClusters}{$superC};
			#print "At cutoff=$Cutoff\tnsubc = $nsubc\t nsubc_previous=$nsubcPre\n";
			#print "\tnew cluster\n\t";
			#print Dumper $Clusters{SubClusters}{$superC};
			
						
			if ($nsubc != $nsubcPre) {
				my @OrderedS = sort { $Clusters{SubClusters}{$superC}{$a} <=> $Clusters{SubClusters}{$superC}{$b} } keys $Clusters{SubClusters}{$superC};
				my @OrderedE = @{$Clusters{SubClusters}{$superC}}{@OrderedS};
				foreach my $start (@OrderedS) {
				push @ClusterArray, [$start,$Clusters{SubClusters}{$superC}{$start},$Cutoff];
				}
				$nsubcPre = $nsubc;
			}
			elsif ($nsubc == $nsubcPre && defined $TempSubClusterRef) {
				# print "\tinside the second loop\n";
				# print "\tTemporary cluster\n\t";
				# print Dumper $TempSubClusterRef;
				my $OverlapCheck = 1;
				foreach my $subC (keys $Clusters{SubClusters}{$superC}) {
					foreach my $subCpre (keys $TempSubClusterRef) {
						my ($ss, $se, $es, $ee);
						$ss = $subC - $subCpre;
						$se = $subC - $TempSubClusterRef->{$subCpre};
						$es = $Clusters{SubClusters}{$superC}{$subC} - $subCpre;
						$ee = $Clusters{SubClusters}{$superC}{$subC} - $TempSubClusterRef->{$subCpre};
						#print "\t\tFor subC=$subC, subCpre=$subCpre: ($ss, $es, $se, $ee)";
						if (($ss*$es)> 0 && ($se*$ee) > 0) {
							$OverlapCheck = 0; # one of new sub clusters does not overlap with the previous sub cluster in consideration, check the next previous sub cluster
							#print "\tOverlapCheck = $OverlapCheck  continue checking\n";
						}
						else {
							$OverlapCheck = 1;
							#print "\tOverlapCheck = $OverlapCheck  stop checking, found one overlap\n";
							last; # one of new sub clusters overlaps with one of the previous sub clusters, stop checking
						}
					}
					if ($OverlapCheck == 0) {
						#print "\tOverlapCheck=0 passed subC=$subC\n";
						last; # one of new sub clusters does not overlap with any of the previous sub clusters, stop checking for other non-overlapping new clusters. One is enough!
					}
				}
				if ($OverlapCheck == 0) {
					my @OrderedS = sort { $Clusters{SubClusters}{$superC}{$a} <=> $Clusters{SubClusters}{$superC}{$b} } keys $Clusters{SubClusters}{$superC};
					my @OrderedE = @{$Clusters{SubClusters}{$superC}}{@OrderedS};
					foreach my $start (@OrderedS) {
					push @ClusterArray, [$start,$Clusters{SubClusters}{$superC}{$start},$Cutoff];
					}
				}
			}
		}
		$TempSubClusterRef = $Clusters{SubClusters}{$superC};
		delete $Clusters{SubClusters}{$superC}; # since we track the number of subclusters by nsubc and nsubcPre we delete the hash--> important to correctly track the number
		#print "\n\n";
	}
	if ($Clusters{SuperClusters}{$superC} - $superC >= $MinPts) {
		push @ClusterArray, [$superC,$Clusters{SuperClusters}{$superC}-1,9.5]; # Artificially recording the super cluster (-1 b/c end has included the start of the next)
	}
	#print Dumper \%Clusters;
}
#print Dumper \@ClusterArray;



########################  Writing to files  ##############################

my $OutFile1 = "./Results/$inputFile.SuperClusters.plot";
open (OUT, ">$OutFile1");

for (my $i = 0; $i < scalar @ClusterArray; $i++) {
	print OUT "$ClusterArray[$i][0]\t$ClusterArray[$i][1]\t$ClusterArray[$i][2]\n";
}
close (OUT);

########################################

my $OutFile2 = "./Results/$inputFile.clusters";

open (OUT, ">$OutFile2");

print OUT "Cluster\tGene/Drug\tMutation/Gene\tDegree_Connectivity\tCloseness_Centrality\tGeodesic_From_Centroid\tRecurrence\n";

for (my $i = 0; $i < scalar @ClusterArray; $i++) {
	for (my $j = $ClusterArray[$i][0]; $j <= $ClusterArray[$i][1]; $j++) {
		my @characters = split(":",$InitialSet[$j][0]);
		my $Gene = $characters[0];
		my $Mutation = $characters[1];
		print OUT "$i\t$Gene\t$Mutation\t0\t0\t1\t0\n";
	}
}

close (OUT);


############################  Plotting  #################################

system ("Rscript BruteForceAutoClustersLines.R $inputFile ./Results/$inputFile.SuperClusters.plot SuperClusters.$inputFile.pdf $Epsilon $MinPts");
print "Done.\n";

#########################################################################