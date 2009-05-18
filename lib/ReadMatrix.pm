package ReadMatrix;
use strict;
use warnings;
use 5.010;

use Exporter qw(import);

our @EXPORT_OK = qw(real_matrix);

sub real_matrix  {
    my $str = shift;
    $str =~ m/\[(\d+),(\d+)\]/g
        or die "Matrix format not recognized (missing dimension marker like [4,4]\n";
    my ($w, $h) = ($1, $2);
    print pos($str), $/;
    my $par = qr{ \(  [^()]+ \) }x;
    $str =~ m/\G\(($par(?:,$par)*)\)/g or die "Can't match nested parenthesis";
    $str = $1;
    my @m = split /[(),]+/, $str;
    # leading item is empty, because the string starts with '('
    shift @m;
    my @matrix;
    while (@m) {
        push @matrix, [splice @m, 0, $w];
    }
    return \@matrix;
}


1;
# vim: sw=4 ts=4 expandtab tw=0
