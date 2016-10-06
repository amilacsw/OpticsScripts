use strict;
use warnings;

use LWP::Simple;
use FileHandle;
use Data::Dumper qw(Dumper);
use Scalar::Util qw(looks_like_number);
use List::Util qw[min max];

my %SetOfNodes;
my $file = "./Test/amila.test.pairwise";
open(IN, "<$file") || die "Can't open $file: $!";

while (my $line = <IN>) {
	chomp $line;
	my @tabs = split(/\t/,$line);
	my @char19 = split("",$tabs[19]);
	my $dis = $char19[0].$char19[1].$char19[2].$char19[3].$char19[4];
	my $key1 = CombineWords($tabs[0],$tabs[4]);
	my $value1 = CombineWords($tabs[9],$tabs[13]);

	$SetOfNodes{$key1}{distances}{$value1} = $dis;
	$SetOfNodes{$value1}{distances}{$key1} = $dis;

	# if (defined $value1){
	# 	print "$value1\n";
	# }
	# else {print "undefined undefined ------------- undefined\n";}
	#print "$value1\n";
}

foreach my $i (keys %SetOfNodes) {
	$SetOfNodes{$i}{processInfo} = "False";
}

print Dumper \%SetOfNodes;

sub CombineWords {
	my ($word1,$word2)=@_;
	return $word1.":".$word2;
}