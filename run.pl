#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;
use Data::Dumper;

my @hosts = glob "wvbh07{0,2,3,3,5,6,7,8} wvbh06{6,9} wthp00{6,9} wthp010 wthp10{4,4,5,5,6,6}";
my $parallel_jobs = @hosts;
my $revoke;
$revoke = 1 if $ARGV[0] && $ARGV[0] eq 'revoke';

my $dir = 'data/' . int(rand() * 10_000);
unless ($revoke) {
    die "bad luck" if -e $dir;
    mkdir $dir or die "Can't mkdir `$dir': $!";

    print "Writing data to `$dir'\n";
}


my $b = 0;

my $pm = Parallel::ForkManager->new($parallel_jobs);

my $count = -1;
for my $phi (0..90) {
    $count++;
    my $pid = $pm->start and next;
    my $host = $hosts[$count % @hosts];
    if ($revoke) {
        system 'ssh', '-x', $host, 'killall', 'cppspin';
    } else {
        my $fn = sprintf "%s/bz%+.2f,phi%02d.dat", $dir, $b, $phi;
        my @args = (
            -b => $b,
            -e => 0.02,
            -o => $fn,
            -r => 0.05,
            -p => $phi,
            '-q',
            -n => 19,
        );
        print "START: ($host) Phi = ", $phi, " degrees\n";
        my $ret = system('ssh', '-x', $host, './run.sh', @args);
        print "END:   ($host) Phi = ", $phi, " degrees\n";
        if ($ret != 0) {
            warn "can't run ssh $host run.sh: $?, $!\n";
            warn "re-running it locally...\n";
            system './cppspin', @args;
        }
    }
    $pm->finish;
}
$pm->wait_all_children;

print "finished run `$dir'\n";

# vim: ft=perl sw=4 ts=4 expandtab
