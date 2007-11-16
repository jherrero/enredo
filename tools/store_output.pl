#!/usr/bin/env perl
use warnings;
use strict;

use Getopt::Long;
use Sys::Hostname;
use Cwd qw(realpath);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::MethodLinkSpeciesSet;
use Bio::EnsEMBL::Compara::SyntenyRegion;
use Bio::EnsEMBL::Compara::DnaFragRegion;

my $description = q'
PROGRAM: store_ouptut.pl

DESCRIPTION: This software allows to store the output of Enredo in an Ensembl
  Compara database

SYNOPSIS:

perl store_output.pl [options] -i enredo.out

OPTIONS:
 --reg_url mysql://user:pass@host:port [default: -none-]
 --compara NAME [default: Multi]
 --db_url mysql://user:pass@host:port/compara_db_name
 --method_link_type TYPE [default: ENREDO]
 --name method_link_species_set.name [default: derived from
      the species names and the method_link_type]
 --source method_link_species_set.source [default: ensembl]
 --url method_link_species_set.url [default: full path to input file]
 --force_strand [default: -disabled-]
 --help [default: -disabled-]

* reg_url: uses a URL to auto-load the registry (needed to connect to the
    core databases)
* compara: name of the database in the Registry (should be Multi or compara)
* db_url: alternatively you can use this way to define the compara
    database you want to use to store the results
* method_link_type: method_link_type for the new data
* name: the name for the new MethodLinkSpeciesSet entry. A name will
    be automatically created from the name of the species and the
    method_link_type if none is provided
* source: the source for the new MethodLinkSpeciesSet entry. Default values
    is "ensembl".
* url: the URL for the new MethodLinkSpeciesSet entry. Default value is
    the location of the input_file.
* force_strand: palindromic segments may have an undefined
    strand. This option sets these strands to 1 (you should probably avoid
    this option and use the EstimateStrand.pl script instead)
* estimate_strand: palindromic segments may have an undefined
    strand. This option requires the EstimateStrand.pl script which
    uses PrePecan to get pariwise alignments in both directions and
    guesses which one is the best
* run_estimate_strand: run string for EstimateStrand.pl script. Default is:
    "perl ./EstimateStrand.pl"
* help: show this help message

EXAMPLE:
  perl store_output.pl --reg_url mysql://user@ens-livemirror \
      --db_url mysql://user:pass@compara1/my_compara_db -i enredo.out
';

my $reg_url;
my $db_url;
my $compara = "Multi";
my $method_link_type = "ENREDO";
my $name = undef;
my $source = "ensembl";
my $url = "";
my $input_file;
my $force_strand = 0;
my $estimate_strand = 0;
my $run_estimate_strand = "perl ./EstimateStrand.pl";
my $dry_run = 0;
my $help = 0;

GetOptions(
    "method_link_type=s" => \$method_link_type,
    "name=s" => \$name,
    "source=s" => \$name,
    "url=s" => \$name,
    "reg_url=s" => \$reg_url,
    "compara=s" => \$compara,
    "db_url=s" => \$db_url,
    "i=s" => \$input_file,
    "force_strand" => \$force_strand,
    "estimate_strand" => \$estimate_strand,
    "run_estimate_strand=s" => \$run_estimate_strand,
    "dry_run|dry-run|n" => \$dry_run,
    "help" => \$help,
  );

if ($help or !$input_file) {
  print $description;
  exit();
}

my $reg = "Bio::EnsEMBL::Registry";
$reg->no_version_check(1);
if ($reg_url) {
  $reg->load_registry_from_url($reg_url);
}
if ($db_url =~ /mysql\:\/\/([^\@]+\@)?([^\:\/]+)(\:\d+)?(\/\w+)?/) {
  my $user_pass = $1;
  my $host = $2;
  my $port = $3;
  my $dbname = $4;

  $user_pass =~ s/\@$//;
  my ($user, $pass) = $user_pass =~ m/([^\:]+)(\:.+)?/;
  $pass =~ s/^\:// if ($pass);
  $port =~ s/^\:// if ($port);
  $dbname =~ s/^\/// if ($dbname);

  Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
        -host=> $host,
        -user => $user,
        -pass => $pass,
        -port => $port,
        -species => $dbname,
        -dbname => $dbname,
    );
  $compara = $dbname;
}

my ($all_synteny_regions, $genome_dbs) = parse_file($input_file);
print "## Enredo parser\n";
print "# Found ", scalar(@$all_synteny_regions), " blocks\n";

my $method_link_species_set_adaptor = $reg->get_adaptor(
    $compara, "compara", "MethodLinkSpeciesSet");
if (!$name) {
  $name = join("-", map {$_->name =~ /(.)\S+\s(.{3})/; $1.".".$2} values %$genome_dbs).
      " ".lc($method_link_type);
}
if (!$url) {
  $url = "file://[".hostname."]".realpath($input_file);
}

my $method_link_species_set = Bio::EnsEMBL::Compara::MethodLinkSpeciesSet->new(
        -method_link_type => $method_link_type,
        -species_set => [values %$genome_dbs],
        -method_link_class => "SyntenyRegion.synteny",
        -name => $name,
        -source => $source,
        -url => $url,
    );

if ($dry_run) {
  print "# Method link type: ", $method_link_species_set->method_link_type, "\n";
  print "# GenomeDBs: ", join (" -- ", map {$_->name} @{$method_link_species_set->species_set}), "\n";
  print "# Method link species set name: ", $method_link_species_set->name, "\n";
  print "# Method link species set source: ", $method_link_species_set->source, "\n";
  print "# Method link species set URL: ", $method_link_species_set->url, "\n";
  print "\n  -- DRY-RUN MODE: NO DATA HAVE BEEN STORED! --\n\n!";
  exit(0);
}

$method_link_species_set = $method_link_species_set_adaptor->store($method_link_species_set);
print "# Method link type: ", $method_link_species_set->method_link_type, "\n";
print "# GenomeDBs: ", join (" -- ", map {$_->name} @{$method_link_species_set->species_set}), "\n";
print "# Method link species set ID: ", $method_link_species_set->dbID, "\n";
print "# Method link species set name: ", $method_link_species_set->name, "\n";
print "# Method link species set source: ", $method_link_species_set->source, "\n";
print "# Method link species set URL: ", $method_link_species_set->url, "\n";

my $synteny_region_adaptor = $reg->get_adaptor(
    $compara, "compara", "SyntenyRegion");
foreach my $this_synteny_region (@$all_synteny_regions) {
  $this_synteny_region->method_link_species_set_id($method_link_species_set->dbID);
  $synteny_region_adaptor->store($this_synteny_region);
}
exit(0);


=head2 parse_file

  Arg[1]        string $filename
  Example       my $blocks = parse_file("enredo.out");
  Description:  reads the input file and returns a ref to an array of
                hashes containing the syntenic blocks as read from the
                enredo output file
  ReturnType:   in array context:
                  ([Bio::EnsEMBL::Compara::SyntenyRegion],
                    {Bio::EnsEMBL::Compara::GenomeDB})
                in scalar context:
                  [Bio::EnsEMBL::Compara::SyntenyRegion]

=cut


sub parse_file {
  my ($file_name) = @_;

  my $synteny_regions = [];
  my $dnafrag_regions = [];
  my $description = undef;

  my $genome_db_adaptor = $reg->get_adaptor($compara, "compara", "GenomeDB");
  my $dnafrag_adaptor = $reg->get_adaptor($compara, "compara", "DnaFrag");
  my $dnafrag_region_adaptor = $reg->get_adaptor($compara, "compara", "DnaFragRegion");
  my $dnafrags;
  my $null_strands = 0;

  open(FILE, $file_name) or return undef;
  while(<FILE>) {
    next if (/^#/);
    if (/^$/) {
      #new block
      if (@$dnafrag_regions) {
        my $this_synteny_region = Bio::EnsEMBL::Compara::SyntenyRegion->new();
        if ($estimate_strand and grep {$_ == 0} map {$_->dnafrag_strand} @$dnafrag_regions) {
          my $fasta_to_dnafrag_region = {};
          my $run_str = $run_estimate_strand;
#           for (my $i = 0; $i < @$dnafrag_regions; $i++) {
#             my $this_dnafrag_region = $dnafrag_regions->[$i];
#             $this_dnafrag_region->adaptor($dnafrag_region_adaptor);
#             print $this_dnafrag_region->dnafrag_strand, " ", $this_dnafrag_region->slice->name, "\n";
#           }
          for (my $i = 0; $i < @$dnafrag_regions; $i++) {
            my $this_dnafrag_region = $dnafrag_regions->[$i];
            $this_dnafrag_region->adaptor($dnafrag_region_adaptor);
            my $filename = "seq_EstStr_".$$."_$i.fa";
            $fasta_to_dnafrag_region->{$filename} = $this_dnafrag_region;
            open(FASTA_FILE, ">$filename") or die;
            print FASTA_FILE ">", $this_dnafrag_region->slice->name, "\n";
            my $seq = $this_dnafrag_region->slice->get_repeatmasked_seq(undef,1)->seq;
            $seq =~ s/(.{60})/$1\n/g;
            print FASTA_FILE $seq, "\n";
            close(FASTA_FILE);
            if ($this_dnafrag_region->dnafrag_strand == 0) {
              $run_str .= " --query $filename";
            } else {
              $run_str .= " --known $filename";
            }
          }
#           print "\n\nRUNNING\n$run_str\n";
#           print "[enter]";
#           <STDIN>;
          my $result = qx"$run_str 2> /dev/null";
          if (!$result) {
            die "====================================================\n".
                "The following command line:\n".
                "  $run_str\n".
                "did not return any output\n".
                "Use the --run_estimate_strand option to change its\n".
                "default value ($run_estimate_strand)\n".
                "====================================================\n";
          }
#           print "\n\nRESULTS\n$result\n";
          foreach my $line (split("\n", $result)) {
#             print "LINE: $line\n";
            my ($filename, $strand, $score) = ($line =~ /(.+) (.+) ([\-\d.e]+)/);
            if (defined($fasta_to_dnafrag_region->{$filename})) {
              if ($strand eq "FWD") {
                $fasta_to_dnafrag_region->{$filename}->dnafrag_strand(1);
              } elsif ($strand eq "RVS") {
                $fasta_to_dnafrag_region->{$filename}->dnafrag_strand(-1);
              } else {
                warn("Error");
              }
            } else {
              warn("Error: $line; $filename\n");
            }
          }
          foreach my $filename (keys %$fasta_to_dnafrag_region) {
            unlink($filename);
          }
          unlink("output.mfa");
#           for (my $i = 0; $i < @$dnafrag_regions; $i++) {
#             my $this_dnafrag_region = $dnafrag_regions->[$i];
#             print $this_dnafrag_region->dnafrag_strand, " ", $this_dnafrag_region->slice->name, "\n";
#           }
#           print "[enter]";
#           <STDIN>;
        }
        foreach my $this_dnafrag_region (@$dnafrag_regions) {
          $this_synteny_region->add_child($this_dnafrag_region);
        }
        $this_synteny_region->{_enredo_description} = $description;
        push(@$synteny_regions, $this_synteny_region);
      }
      $dnafrag_regions = [];
      $description = undef;
    } elsif (/^block( .*)?$/) {
      if (@$dnafrag_regions) {
        die "Starting a new block before ending previous one...\n";
      }
      $description = $1;
    } elsif (/^([^\:]+):([^\:]+):(\d+):(\d+) \[(.*)\]/) {
      my $species = $1;
      my $chromosome = $2;
      my $start = $3;
      my $end = $4;
      my $strand = $5;
      if (!$genome_dbs->{$species}) {
        $genome_dbs->{$species} = $genome_db_adaptor->fetch_by_registry_name($species);
        die if (!$genome_dbs->{$species});
      }
      if (!$dnafrags->{$species}->{$chromosome}) {
        $dnafrags->{$species}->{$chromosome} = $dnafrag_adaptor->
            fetch_by_GenomeDB_and_name($genome_dbs->{$species}, $chromosome);
        die if (!$dnafrags->{$species}->{$chromosome});
      }
      if ($strand == 0) {
        $null_strands++;
        if ($force_strand) {
          $strand = 1;
        }
      }
      my $this_dnafrag_region = Bio::EnsEMBL::Compara::DnaFragRegion->new(
          -dnafrag_id => $dnafrags->{$species}->{$chromosome}->dbID,
          -dnafrag_start => $start,
          -dnafrag_end => $end,
          );
      $this_dnafrag_region->dnafrag_strand($strand);
      push(@$dnafrag_regions, $this_dnafrag_region);
    } else {
      die "Cannot understand line: $_";
    }
  }
  if (@$dnafrag_regions) {
    my $this_synteny_region = Bio::EnsEMBL::Compara::SyntenyRegion->new();
    foreach my $this_dnafrag_region (@$dnafrag_regions) {
      $this_synteny_region->add_child($this_dnafrag_region);
    }
    $this_synteny_region->{_enredo_description} = $description;
    push(@$synteny_regions, $this_synteny_region);
  }
  close(FILE);

  ## Warning about null strands
  if ($null_strands) {
    if ($estimate_strand) {
      warn("The strand for $null_strands regions has been estimated with PrePecan aligner");
    } elsif ($force_strand) {
      warn("The strand for $null_strands regions has been set to 1");
    } else {
      warn("$null_strands regions have been ignored because the strand was undefined\n".
          "Use the --force_strand or --estimate_strand (recommended) to solve this");
    }
  }

  if (wantarray) {
    return ($synteny_regions, $genome_dbs);
  } else {
    return $synteny_regions;
  }
}

