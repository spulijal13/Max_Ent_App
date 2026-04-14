use strict;
use warnings;

my $inputfile  = $ARGV[0] or die "Usage: perl score3.pl <inputfile> <outputfile>\n";
my $outputfile = $ARGV[1] or die "Usage: perl score3.pl <inputfile> <outputfile>\n";

my $usemaxent = 1;

my %me2x5 = makescorematrix('splicemodels/me2x5');
my %seq   = makesequencematrix('splicemodels/splice5sequences');

open(my $FILE, "<", $inputfile)  or die "can't open input file $inputfile\n";
open(my $OUT,  ">", $outputfile) or die "can't create output file $outputfile\n";

while (<$FILE>) {
    chomp;
    s/\cM//g;

    next if /^\s*$/;
    next if /^>/;

    my ($id, $str) = split(/\s+/, $_, 2);
    next unless defined $id and defined $str;

    $str =~ s/^\s+|\s+$//g;
    $str = uc($str);

    if ($str =~ /[^ACGT]/) {
        print $OUT $id . "\t" . $str . "\tNA\n";
        next;
    }

    if (length($str) != 9) {
        print $OUT $id . "\t" . $str . "\tNA\n";
        next;
    }

    if ($usemaxent) {
        my $rest = getrest($str);

        if (!exists $seq{$rest} || !exists $me2x5{$seq{$rest}}) {
            print $OUT $id . "\t" . $str . "\tNA\n";
            next;
        }

        my $score = sprintf(
            "%.2f",
            log2(scoreconsensus($str) * $me2x5{$seq{$rest}})
        );

        print $OUT $id . "\t" . $str . "\t" . $score . "\n";
    }
    else {
        print $OUT $id . "\t" . $str . "\tNA\n";
    }
}

close($FILE);
close($OUT);

sub makesequencematrix {
    my $file = shift;
    my %matrix;
    my $n = 0;

    open(my $SCOREF, "<", $file) or die "Can't open $file\n";
    while (<$SCOREF>) {
        chomp;
        s/\s//g;
        $matrix{$_} = $n;
        $n++;
    }
    close($SCOREF);

    return %matrix;
}

sub makescorematrix {
    my $file = shift;
    my %matrix;
    my $n = 0;

    open(my $SCOREF, "<", $file) or die "Can't open $file\n";
    while (<$SCOREF>) {
        chomp;
        s/\s//g;
        $matrix{$n} = $_;
        $n++;
    }
    close($SCOREF);

    return %matrix;
}

sub getrest {
    my $seq = shift;
    my @seqa = split(//, uc($seq));
    return $seqa[0] . $seqa[1] . $seqa[2] . $seqa[5] . $seqa[6] . $seqa[7] . $seqa[8];
}

sub scoreconsensus {
    my $seq = shift;
    my @seqa = split(//, uc($seq));

    my %bgd = (
        'A' => 0.27,
        'C' => 0.23,
        'G' => 0.23,
        'T' => 0.27,
    );

    my %cons1 = (
        'A' => 0.0040,
        'C' => 0.0032,
        'G' => 0.9896,
        'T' => 0.0032,
    );

    my %cons2 = (
        'A' => 0.0034,
        'C' => 0.0039,
        'G' => 0.0042,
        'T' => 0.9884,
    );

    my $addscore = $cons1{$seqa[3]} * $cons2{$seqa[4]} /
                   ($bgd{$seqa[3]} * $bgd{$seqa[4]});

    return $addscore;
}

sub log2 {
    my ($val) = @_;
    return log($val) / log(2);
}