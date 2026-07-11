#!/usr/bin/perl
use strict;
use warnings;

sub usage {
    my ($handle) = @_;
    print {$handle} <<'EOF';
Usage: localMaxSkmelTADN.pl INPUT.tsv

Select one row per consecutive gene group from a tab-separated promoter table.
Rows are ranked by the sum of the SK-Mel-29 DN and TA p73 activity columns;
ties prefer matching motif/gene direction. Input must contain at least 24
columns and be grouped by gene. Selected original rows are written to stdout.

Options:
  -h, --help  Show this help and exit
EOF
}

if (@ARGV == 1 && ($ARGV[0] eq '-h' || $ARGV[0] eq '--help')) {
    usage(*STDOUT);
    exit 0;
}

if (@ARGV != 1) {
    usage(*STDERR);
    exit 2;
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
