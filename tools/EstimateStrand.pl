#!/usr/bin/env perl

use warnings;
use strict;

use Bio::SeqIO;
use Bio::AlignIO;
use Getopt::Long;

my $program_name = "estimate_strand.pl";
my $version = "0.0.1a";

my $defaults = {
  java => "java",
  pecan_jar => "/home/jherrero/Downloads/pecan.0.6.jar",
  prepecan_prefix => "",
};

my $JAVA = undef;
my $PECAN_JAR = undef;
my $PREPECAN_PREFIX = undef;

my $description = qq!
  SYNOPSIS:
    perl $program_name [configuration_options] \
        --known known_seq.fa --query query_seq.fa
  WHERE
    stranded_seq.fa is a sequence in fasta format. This program will
    tell the user whether query_seq is better aligned to stranded_seq.fa
    on the forward or the reverse strand. You can have several known and
    several query sequences. Each query sequence will be tested against
    all the known sequence in order to deciede about the strand. If no
    known sequences are defined and there are several query sequences,
    the first one will be used in the forward (current) strand and all
    subsequent ones will be stranded according to this one.
  CONFIGURATION:
    --java PATH_TO_JAVA [default: !.$defaults->{java}.qq!]
    --pecan_jar PATH_TO_PECAN_JAR [default: !.$defaults->{pecan_jar}.qq!]
    --prepecan_prefix PREPECAN_PREFIX [default: JAVA -server -cp PECAN_JAR bp.pecan.utils.PrePecan]
!;

sub help {
  print "$program_name $version\n";
  print $description;
}

my $help;
my $known_seqs = [];
my $query_seqs = [];

GetOptions(
    "help" => \$help,
    "known=s" => \@$known_seqs,
    "query=s" => \@$query_seqs,
    "java=s" => \$JAVA,
    "pecan_jar=s" => \$PECAN_JAR,
    "prepecan_prefix=s" => \$PREPECAN_PREFIX,
  );

$JAVA = $defaults->{java} if (!$JAVA);
$PECAN_JAR = $defaults->{pecan_jar} if (!$PECAN_JAR);
$PREPECAN_PREFIX = "$JAVA -server -cp $PECAN_JAR bp.pecan.utils.PrePecan" if (!$PREPECAN_PREFIX);

if ($help) {
  help();
  exit(0);
}

if (!@$query_seqs) {
  print STDERR "No query sequences. Please look for help: $program_name --help\n";
} elsif (!@$known_seqs and @$query_seqs == 1) {
  print STDERR "No known sequences and only 1 query sequence. Please look for help: $program_name --help\n";
}

if (!@$known_seqs) {
  my $this_query_seq = shift(@$query_seqs);
  print "$this_query_seq FWD 100000\n";
  push(@$known_seqs, $this_query_seq);
}

foreach my $this_query_seq (@$query_seqs) {
#   print "QUERY: $this_query_seq\n";
  my $query_in = new Bio::SeqIO(-file => $this_query_seq, -format => 'fasta');
  my $seq = $query_in->next_seq();
  my $rev = $seq->revcom;
  my $query_out = new Bio::SeqIO(-file => ">".$this_query_seq."_revcom", -format => 'fasta');
  $query_out->write_seq($rev);
  my $score = 0;
  foreach my $this_known_seq (@$known_seqs) {
#     print "  KNOWN: $this_known_seq\n";
    my $identity_fwd = get_percentage_identity($this_known_seq, $this_query_seq);
    $score += $identity_fwd;

    my $identity_rvs = get_percentage_identity($this_known_seq, $this_query_seq."_revcom");
    $score -= $identity_rvs;
  }
  if ($score >= 0) {
    print "$this_query_seq FWD $score\n";
  } else {
    print "$this_query_seq RVS $score\n";
  }
  unlink($this_query_seq."_revcom");
}

sub get_percentage_identity {
  my ($first_file, $second_file) = @_;

  my $run_str = join(" ", $PREPECAN_PREFIX, "-F", $first_file, $second_file, "-E", "'(1,2);'", "&>", "/dev/null");
  print STDERR "$run_str\n";
  my $result = qx"$run_str";
  if (!-e "output.mfa") {
    die "Error while running PrePecan: output.mfa not found!\n";
  }
  my $percentage_identity = evaluate_alignment("output.mfa");

  return $percentage_identity;
}

sub evaluate_alignment {
  my ($alignment_file, $format) = @_;

  $alignment_file ||= "output.mfa";
  $format ||= "fasta";
  my $pecan_in = Bio::AlignIO->new(-file => $alignment_file, -format => $format);
  my $aln = $pecan_in->next_aln(100); # There sohould be one single alignment only!
  my $percentage_identity = $aln->overall_percentage_identity;

  return $percentage_identity;
}

# 1979
# human:10:5187882:5297899 [-1] l=110018
# rat:17:77153557:77205935 [0] l=52379

# 25055
# rat:1:162333305:162452957 [-1] l=119653
# human:11:6085928:6106795 [0] l=20868
# mouse:7:105149660:105162313 [0] l=12654
# mouse:7:105162259:105188221 [0] l=25963
# dog:21:32647610:32666600 [0] l=18991

# 26179
# rat:6:105230387:105285333 [0] l=54947
# dog:8:47145746:47157897 [0] l=12152
# mouse:12:82297738:82340252 [1] l=42515
