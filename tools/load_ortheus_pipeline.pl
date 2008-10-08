#!/usr/bin/env perl

use warnings;
use strict;

=head1 NAME

load_ortheus_pipeline.pl

=head1 AUTHORS

 Javier Herrero (jherrero@ebi.ac.uk)

=head1 COPYRIGHT

This script is part of the Ensembl project http://www.ensembl.org

=head1 DESCRIPTION

This script is intended to load all the pieces that are neede to start
an Ortheus pipeline out of a set of SyntenyRegion entries, ie, the output
of Enredo.

=head1 SYNOPSIS

perl load_ortheus_pipeline.pl --help

perl load_ortheus_pipeline.pl
    <db_url=s>
    <species_tree|tree=s>
    [mlssid|mlss_id|method_link_species_set_id=s]
    [core_db=s]

=head1 REQUIREMENTS

This script uses mysql and the Ensembl Perl API.

=head1 ARGUMENTS

=head2 GETTING HELP

=over

=item B<[--help]>

  Prints help message and exits.

=back

=head2 PARAMETERS

=over

=item B<--db_url mysql://user[:passwd]@host[:port]/dbname>

Mandatory

URL for the compara database containing the synteny_region data. The pipeline will
be run on this instance. Also, this script expects this database to contain all the
tables needed for running the Hive system:

 ensembl-compara/sql/table.sql
 ensembl-compara/sql/pipeline-tables.sql
 ensembl-hive/sql/tables.sql


=item B<--species_tree "((22:0.1381,(25:0.0770,3:0.0817):0.2526):0.0230,(39:0.1477,45:0.1592):0.0394);">

Mandatory

A tree in Newick format where the name of the leaves are the genome_db_ids of the species to be aligned.

=item B[--mlss_id 50001]

Optional

The method_link_species_set_id for the SyntenyRegions to be aligned. If none is specified, the script
will propose the user to pick uk one from all the available ones that are relevant.

=item B[--core_db my_db_name]

Optional

The name of the database for storing the ancestral sequences. If none is specified, the script will use
use one based on the name of the compara database: [compara_db_name]_ancestral_core

=back

=head1 INTERNAL METHODS

=cut


use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Compara::GenomeDB;
use Bio::EnsEMBL::Compara::MethodLinkSpeciesSet;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Getopt::Long;

Bio::EnsEMBL::Registry->no_version_check(1);

my $help;
my $db_url;
my $core_db;
my $mlss_id;
my $species_tree;

GetOptions(
    "help" => \$help,
    "db_url=s" => \$db_url,
    "core_db=s" => \$core_db,
    "mlssid|mlss_id|method_link_species_set_id=s" => \$mlss_id,
    "species_tree|tree=s" => \$species_tree,
  );

if ($help or !$db_url or !$species_tree) {
  exec("/usr/bin/env perldoc $0");
}

my $compara_dba = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(-url => $db_url);

my $method_link_species_set = get_method_link_species_set($compara_dba, $mlss_id);
$mlss_id = $method_link_species_set->dbID;

create_tables($compara_dba);

my $ancestral_genome_db = create_core_genome_db($compara_dba);

create_core_database($compara_dba, $ancestral_genome_db, $core_db);

my $analysis = create_analysis($compara_dba);

my $ortheus_mlss = create_ortheus_mlss($compara_dba, $method_link_species_set);

create_analysis_jobs($compara_dba, $method_link_species_set, $analysis, $ortheus_mlss, $species_tree);

sub get_method_link_species_set {
  my ($compara_dba, $mlss_id) = @_;
  my $method_link_species_set;

  if ($mlss_id) {
    $method_link_species_set = $compara_dba->get_MethodLinkSpeciesSetAdaptor->fetch_by_dbID($mlss_id);
  }

  while (!$method_link_species_set) {
    print "Here is the list of all relevant MethodLinkSpeciesSet in the database:\n";
    my @method_link_species_sets = grep {!$_->method_link_class or $_->method_link_class =~ /^SyntenyRegion/}
        @{$compara_dba->get_MethodLinkSpeciesSetAdaptor->fetch_all()};
    foreach my $this_mlss (@method_link_species_sets) {
      next if ($this_mlss->method_link_class and $this_mlss->method_link_class !~ /^SyntenyRegion/);
      print $this_mlss->dbID, ". ", $this_mlss->name, ":\n    ", $this_mlss->method_link_type, " - ",
          join(", ", map {$_->name} @{$this_mlss->species_set}), "\n";
    }
    print "Please, enter the ID for the MethodLinkSpeciesSet you want to use: ";
    $mlss_id = <STDIN>;
    $method_link_species_set = $compara_dba->get_MethodLinkSpeciesSetAdaptor->fetch_by_dbID($mlss_id);
  }

  return $method_link_species_set;
}


sub create_tables {
  my ($compara_dba) = @_;

  my $create_genomic_align_tree_table = qq!
      CREATE TABLE IF NOT EXISTS `genomic_align_tree` (
        `node_id` bigint(20) unsigned NOT NULL auto_increment,
        `parent_id` bigint(20) unsigned NOT NULL default '0',
        `root_id` bigint(20) unsigned NOT NULL default '0',
        `left_index` int(10) NOT NULL default '0',
        `right_index` int(10) NOT NULL default '0',
        `distance_to_parent` double NOT NULL default '1',
        PRIMARY KEY  (`node_id`),
        KEY `parent_id` (`parent_id`),
        KEY `root_id` (`root_id`),
        KEY `left_index` (`left_index`),
        KEY `right_index` (`right_index`)
      ) ENGINE=MyISAM CHARSET=latin1
  !;

  $compara_dba->dbc->do($create_genomic_align_tree_table);
}


sub create_core_genome_db {
  my ($compara_dba) = @_;

  my $insert_core_genome_db = qq!
      INSERT IGNORE INTO genome_db (name, taxon_id, assembly)
      VALUES ("Ancestral sequences", 0, "")
  !;

  $compara_dba->dbc->do($insert_core_genome_db);

  return $compara_dba->get_GenomeDBAdaptor->fetch_by_name_assembly("Ancestral sequences");
}


sub create_core_database {
  my ($compara_dba, $ancestral_genome_db, $core_db) = @_;

  if (!$core_db) {
    $core_db = $compara_dba->dbc->dbname."_ancestral_core";
  }

  my $value = $compara_dba->dbc->db_handle->selectrow_array("SHOW DATABASES LIKE \"$core_db\"");
  if (!$value) {
    my $create_core_database = "CREATE DATABASE IF NOT EXISTS $core_db";

    $compara_dba->dbc->do($create_core_database);

    my $table_sql_file;
    foreach my $dir (@INC) {
      if (-e "$dir/Bio/EnsEMBL/Registry.pm" and -e "$dir/../sql/table.sql") {
        $table_sql_file = "$dir/../sql/table.sql";
        last;
      }
    }
    my $exec = "mysql";
    $exec .= " -u".$compara_dba->dbc->username if ($compara_dba->dbc->username);
    $exec .= " -p".$compara_dba->dbc->password if ($compara_dba->dbc->password);
    $exec .= " -h".$compara_dba->dbc->host if ($compara_dba->dbc->host);
    $exec .= " -P".$compara_dba->dbc->port if ($compara_dba->dbc->port);
    $exec .= " $core_db < $table_sql_file";
    system($exec);
  }

  my $locator = "Bio::EnsEMBL::DBSQL::DBAdaptor/host=".
      $compara_dba->dbc->host.
      ";port=".$compara_dba->dbc->port.
      ";user=".$compara_dba->dbc->username.
      ";pass=".$compara_dba->dbc->password.
      ";dbname=".$core_db.
      ";species=Ancestral sequences;disconnect_when_inactive=1";
  my $update_core_genome_db = qq!
      UPDATE genome_db SET locator = "$locator"
      WHERE name = "Ancestral sequences"
  !;

  $compara_dba->dbc->do($update_core_genome_db);

  return $core_db;
}


sub create_analysis {
  my ($compara_dba) = @_;

  my $analysis = new Bio::EnsEMBL::Analysis(
      -logic_name => "Ortheus",
      -module => "Bio::EnsEMBL::Compara::Production::GenomicAlignBlock::Ortheus",
      -parameters => "{max_block_size=>1000000,java_options=>'-server -Xmx1000M',}",
    );
  $compara_dba->get_AnalysisAdaptor()->store($analysis);

  return $analysis;
}


sub create_ortheus_mlss {
  my ($compara_dba, $method_link_species_set) = @_;

  my $ortheus_mlss = new Bio::EnsEMBL::Compara::MethodLinkSpeciesSet(
      -method_link_type => "ORTHEUS",
      -species_set => $method_link_species_set->species_set
    );

  $compara_dba->get_MethodLinkSpeciesSet->store($ortheus_mlss);

  return $ortheus_mlss;
}


sub create_analysis_jobs {
  my ($compara_dba, $method_link_species_set, $analysis, $ortheus_mlss, $species_tree) = @_;

  my $insert_species_tree = qq!INSERT INTO analysis_data(data) VALUES("$species_tree")!;
  my $sth = $compara_dba->dbc->prepare($insert_species_tree);
  $sth->execute;
  my $tree_analysis_data_id = $sth->{'mysql_insertid'};
  my $analysis_id = $analysis->dbID;
  my $ortheus_mlss_id = $ortheus_mlss->dbID;

  my $insert_new_jos = qq!
    INSERT IGNORE INTO analysis_job(analysis_id, input_id)
    SELECT
      $analysis_id,
      CONCAT("{synteny_region_id=>",synteny_region_id,
        ",method_link_species_set_id=>$ortheus_mlss_id,tree_analysis_data_id=>'$tree_analysis_data_id',}")
    FROM synteny_region
    WHERE method_link_species_set_id = !. $method_link_species_set->dbID;

  $compara_dba->dbc->do($insert_new_jos);
}

