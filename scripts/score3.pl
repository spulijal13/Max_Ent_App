use strict;
use warnings;

my $inputfile  = $ARGV[0] or die "Usage: perl score3.pl <inputfile> <outputfile>\n";
my $outputfile = $ARGV[1] or die "Usage: perl score3.pl <inputfile> <outputfile>\n";

my $usemaxent = 1;

my @metables = makemaxentscores();

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

    if (length($str) != 23) {
        print $OUT $id . "\t" . $str . "\tNA\n";
        next;
    }

    if ($usemaxent) {
        my $score = sprintf(
            "%.2f",
            log2(scoreconsensus($str) * maxentscore(getrest($str), \@metables))
        );

        print $OUT $id . "\t" . $str . "\t" . $score . "\n";
    }
    else {
        print $OUT $id . "\t" . $str . "\tNA\n";
    }
}

close($FILE);
close($OUT);

sub hashseq {
    my $seq = shift;
    $seq = uc($seq);
    $seq =~ tr/ACGT/0123/;
    my @seqa = split(//, $seq);
    my $sum = 0;
    my $len = length($seq);
    my @four = (1, 4, 16, 64, 256, 1024, 4096, 16384);
    my $i = 0;

    while ($i < $len) {
        $sum += $seqa[$i] * $four[$len - $i - 1];
        $i++;
    }

    return $sum;
}

sub makemaxentscores {
    my $dir = "splicemodels/";
    my @list = (
        'me2x3acc1', 'me2x3acc2', 'me2x3acc3',
        'me2x3acc4', 'me2x3acc5', 'me2x3acc6',
        'me2x3acc7', 'me2x3acc8', 'me2x3acc9'
    );

    my @metables;
    my $num = 0;

    foreach my $file (@list) {
        my $n = 0;
        open(my $SCOREF, "<", $dir . $file) or die "Can't open $dir$file\n";

        while (<$SCOREF>) {
            chomp;
            s/\s//g;
            $metables[$num]{$n} = $_;
            $n++;
        }

        close($SCOREF);
        $num++;
    }

    return @metables;
}

sub maxentscore {
    my $seq = shift;
    my $table_ref = shift;
    my @metables = @$table_ref;
    my @sc;

    $sc[0] = $metables[0]{hashseq(substr($seq,  0, 7))};
    $sc[1] = $metables[1]{hashseq(substr($seq,  7, 7))};
    $sc[2] = $metables[2]{hashseq(substr($seq, 14, 7))};
    $sc[3] = $metables[3]{hashseq(substr($seq,  4, 7))};
    $sc[4] = $metables[4]{hashseq(substr($seq, 11, 7))};
    $sc[5] = $metables[5]{hashseq(substr($seq,  4, 3))};
    $sc[6] = $metables[6]{hashseq(substr($seq,  7, 4))};
    $sc[7] = $metables[7]{hashseq(substr($seq, 11, 3))};
    $sc[8] = $metables[8]{hashseq(substr($seq, 14, 4))};

    my $finalscore = $sc[0] * $sc[1] * $sc[2] * $sc[3] * $sc[4] /
                     ($sc[5] * $sc[6] * $sc[7] * $sc[8]);

    return $finalscore;
}

sub getrest {
    my $seq = shift;
    my $seq_noconsensus = substr($seq, 0, 18) . substr($seq, 20, 3);
    return $seq_noconsensus;
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
        'A' => 0.9903,
        'C' => 0.0032,
        'G' => 0.0034,
        'T' => 0.0030,
    );

    my %cons2 = (
        'A' => 0.0027,
        'C' => 0.0037,
        'G' => 0.9905,
        'T' => 0.0030,
    );

    my $addscore = $cons1{$seqa[18]} * $cons2{$seqa[19]} /
                   ($bgd{$seqa[18]} * $bgd{$seqa[19]});

    return $addscore;
}

sub log2 {
    my ($val) = @_;
    return log($val) / log(2);
}