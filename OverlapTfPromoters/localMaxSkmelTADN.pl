#!/usr/bin/perl
use strict;
use warnings;

# Check for input file
if (@ARGV != 1) {
    die "Usage: $0 <input.csv>\n";
}

my $input_file = $ARGV[0];

# Open the input file
open my $fh, '<', $input_file or die "Cannot open file $input_file: $!\n";

# Hash to store the maximum activity line for each gene
my %max_activity;

my $localMax="";
my $localMaxScore = -1;
my $prevGene="";

while (<>) {
    chomp;
    my $line = $_;
    my @columns = split /\t/;  # Assuming CSV is tab-separated

    # Ensure the line has enough columns
    if ( @columns < 23 ) {
        warn "Skipping line due to insufficient columns: $line\n";
        exit 1;
    }

    my $p73_skmel_DN = $columns[15];  # Column 16 (0-based index is 15)
    my $p73_skmel_TA = $columns[17];  # Column 18 (0-based index is 17)

    # Calculate the total activity
    my $tp73_activity = $p73_skmel_DN + $p73_skmel_TA;
    my $directionEqual = 0;
    $directionEqual=1 if $columns[5] eq $columns[23];
    my $gene = $columns[21];  # Column 22 (0-based index is 21)
    my $total_activity = $tp73_activity + $directionEqual;

    #print "Debug: $gene\t$tp73_activity + directionEqual:$directionEqual ($columns[5] == $columns[23])\n";

    if ($gene eq $prevGene) {
        if ($localMaxScore < 0 or $localMaxScore < $tp73_activity) {
            $localMax = $line;
            $localMaxScore = $total_activity;
        } elsif ($localMaxScore == $tp73_activity && $directionEqual>0) {
	    #print STDERR "Happened! ($localMaxScore)\n";
            $localMax = $line;
            $localMaxScore = $total_activity;
	}
    } else {
        print "$localMax\n" if $localMax ne "";
        $localMaxScore = $total_activity;
        $localMax = $line;
        $prevGene = $gene;
    }

}
