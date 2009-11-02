#!/usr/bin/perl
use strict;
use warnings;
use lib 'lib';
use ReadMatrix qw(real_matrix);
use Data::Dumper;
use 5.010;

my @x;
my @d;
my $isfirst = 1;
for my $d (@ARGV) {
    my @files = glob "data/$d/*.dat" or die "No files in `data/$d/*'";

    my $i = 0;
    for my $fn (@files) {
        $fn =~ /alpha([\d.]+)\.dat/ or die "can't work with file name `$fn'";
        my $alpha = $1;
        my $m = tpq_for_file($fn);
        if ($isfirst) {
            push @x, $alpha;
        } else {
            if (abs($x[$i]-$alpha) > 1e5) {
                die "Can't merge datasets with different X axes\n"
                    . "  for the ${i}th set I expected $x[$i] and got $alpha";
            }
        }
        push @{$d[$i]}, $m->[0][2] + $m->[1][2];
    } continue {
        $i++;
    }
    $isfirst = 0;
}
use Data::Dumper;
#print Dumper \@x;
#print Dumper \@d;

for my $i (0..$#x) {
    say $x[$i], ' ', join(' ', @{$d[$i]});

}


sub tpq_for_file {
    my $fn = shift;
    open my $handle, '<', $fn or die "Can't open file `$fn' for reading: $!";
    my $m;
    while (<$handle>) {
        if (/^final tpq/) {
            $m = real_matrix($_);
            last;
        }
    }
    close $handle or warn $!;
    return $m;
}

# vim: sw=4 ts=4 expandtab
