#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV != 1) {
    die "Usage: $0 <input.tsv>\n";
}

my $input_file = $ARGV[0];
open my $fh, '<', $input_file or die "Cannot open file $input_file: $!\n";

my $local_max_line;
my $local_max_activity;
my $local_max_direction_equal = 0;
my $previous_gene;

while (my $line = <$fh>) {
    chomp $line;
    my @columns = split /\t/, $line, -1;

    if (@columns < 24) {
        die "Input line has fewer than 24 tab-separated columns: $line\n";
    }

    my $p73_skmel_dn = $columns[15];
    my $p73_skmel_ta = $columns[17];
    my $tp73_activity = $p73_skmel_dn + $p73_skmel_ta;
    my $direction_equal = $columns[5] eq $columns[23] ? 1 : 0;
    my $gene = $columns[21];

    if (!defined $previous_gene || $gene ne $previous_gene) {
        print "$local_max_line\n" if defined $local_max_line;
        $previous_gene = $gene;
        $local_max_line = $line;
        $local_max_activity = $tp73_activity;
        $local_max_direction_equal = $direction_equal;
    } elsif (
        $tp73_activity > $local_max_activity
        || ($tp73_activity == $local_max_activity
            && $direction_equal > $local_max_direction_equal)
    ) {
        $local_max_line = $line;
        $local_max_activity = $tp73_activity;
        $local_max_direction_equal = $direction_equal;
    }
}

print "$local_max_line\n" if defined $local_max_line;

close $fh or die "Cannot close file $input_file: $!\n";
