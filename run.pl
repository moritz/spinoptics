#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;
use Data::Dumper;

my @hosts = glob "wvbh07{0,1,2,3,3,5,7,8} wvbh06{6,9} wthp009 wthp010 wthp10{4,4,5,5,6,6}";
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
for (my $alpha = 0; $alpha <= 1.5; $alpha += 0.002) {
    $count++;
    my $pid = $pm->start and next;
    my $host = $hosts[$count % @hosts];
    if ($revoke) {
        system 'ssh', '-x', $host, 'killall', 'cppspin';
    } else {
        my $fn = sprintf "%s/alpha%.3f.dat", $dir, $alpha;
        my @args = (
            -b => $b,
            -e => 2.0,
            -o => $fn,
            -r => $alpha,
            -p => 40,
            -n => 19,
            '-q',
        );
        print "START: ($host) alpha = $alpha\n";
        my $ret = system('ssh', '-x', $host, './run.sh', @args);
        print "END:   ($host) alpha = $alpha\n";
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
