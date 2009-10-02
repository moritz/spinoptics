#!/usr/bin/perl
use strict;
use warnings;

my $dir = 'data/' . int(rand() * 10_000);
die "bad luck" if -e $dir;

mkdir $dir or die "Can't mkdir `$dir': $!";

print "Writing data to `$dir'\n";

my $b = 0;

for my $phi (0..90) {
    my $angle = $phi / 180 * 3.14159;
    print "Phi = ", $phi, " degrees\n";
        my $fn = sprintf "%s/bz%+.2f,phi%02d.dat", $dir, $b, $phi;
        system('./cppspin',
                -b => $b,
                -e => -2.5,
                -o => $fn,
                -r => -0.5,
                -p => $angle,
                '-q',
                -n => 19,
            ) == 0
            or die "can't run ./cppsin: $?";
}

print "finished run `$dir'\n";

# vim: ft=perl sw=4 ts=4 expandtab
