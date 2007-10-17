use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Getopt::Long;

=head1 NAME

  assess_graph.pl

=head1 AUTHORS

  Javier Herrero

=head1 CONTACT

  Ensembl developers mailing list at <ensembl-dev@ebi.ac.uk>

=head1 DESCRIPTION

This script is intended to assess the results of Enredo, Mercator or other
software which defines syntenic blocks or co-linear regions. The data must
live in an Ensembl compara database and the assembly, repeat masking and gene
information for each species should be available in Ensembl core databases.
This script performs multiple type of analyses, including:
 - overall coverage
 - N50-like score of the regions (weighted average of the length of the segments)
 - duplication stats
 - gene (protein coding) coverage
 - pseudogenes coverage
 - homologue assessment (not implemented yet)
 - ancient repeats assessment (won't be implemented in the near future)

=head2 overall coverage

This will calculate the total coverage of the co-linear regions. If a particular
region is usd in two co-linear blocks, it will still be counted once.

=head2 n50 scores

This calculates the weighted average length of the blocks for each species. This
is done by getting the length of all the blocks, sorting them, adding up the
lengths until the sum reaches 50% of the total length of all the blocks. The N50
is defined as the length of the last segment which was added in order to reach
the 50% threshold.

=head2 duplication stats

This will look for co-linear blocks that contain more than one region per species.

=head1 SYNOPSIS

  perl assess_graph.pl

=cut

my $usage = qq{
perl assess_graph.pl
  Getting help:
    [--help]

  For the data:
    [--registry mysql://user[:passwd]\@host[:port]/[version]
    --db_url|compara_url mysql://user[:passwd]\@host[:port]/db_name
    --mlss_id|method_link_species_set_id method_link_species_set_id
    [--species species]
        Restict the analysis to this species (must be the binomial name of
        the species)

  For the output:
    [--quick] Skip slow analyses. Will only calculate quick coverage
        (without considering repeat information) and will skip the
        gene coverage stats
    [--output_file filename] stores the output in this file
};

my $help;
my $registry;
my $db_url;
my $method_link_species_set_id;

my $species;
my $output_file = undef;
my $quick = 0;

Bio::EnsEMBL::Registry->no_version_check(1);

GetOptions(
    "help" => \$help,
    "registry=s" => \$registry,
    "db_url|compara_url=s" => \$db_url,
    "mlss_id|method_link_species_set_id=i" => \$method_link_species_set_id,
    "species=s" => \$species,
    "output_file=s" => \$output_file,
    "quick" => \$quick,
  );

# Print Help and exit
if ($help or !$db_url or !$method_link_species_set_id) {
    print $usage;
    exit(0);
}

if ($registry and $registry =~ /^mysql:\/\//) {
  Bio::EnsEMBL::Registry->load_registry_from_url($registry);
} else {
  Bio::EnsEMBL::Registry->load_all($registry);
}

my $compara_dba = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(-url => $db_url);

if ($output_file) {
    open(STDOUT, ">$output_file") or die("Cannot open $output_file");
}

## Get main data
my $method_link_species_set = get_MethodLinkSpeciesSet($compara_dba, $method_link_species_set_id);
my $species_tree = get_tree_by_GenomeDBs($method_link_species_set->species_set);
my $all_synteny_regions = fetch_all_SyntenyRegions($compara_dba, $method_link_species_set);
my $all_dnafrag_regions = get_all_DnaFragRegions($all_synteny_regions);

print_genomic_coverage($all_dnafrag_regions, $quick);
print_N50_stats($all_dnafrag_regions);
print_duplication_stats($all_synteny_regions, $species_tree);
if (!$quick) {
  print_gene_coverage($all_dnafrag_regions, $method_link_species_set, "protein_coding", $species);
  print_gene_coverage($all_dnafrag_regions, $method_link_species_set, "pseudogene", $species);
}

exit(0);


sub print_gene_coverage {
  my ($all_dnafrag_regions, $method_link_species_set, $biotype, $species) = @_;

  my $genome_dbs = $method_link_species_set->species_set;

  foreach my $this_genome_db (@$genome_dbs) {
    next if ($species and $this_genome_db->name ne $species);
    my $gene_adaptor = $this_genome_db->db_adaptor->get_GeneAdaptor();
    my $all_genes = $gene_adaptor->fetch_all_by_biotype($biotype);
    my $sizes = [];
    my $maxs = [];
    my $totals = [];
    my ($fully, $broken, $partially, $uncovered) = (0,0,0);
    print " + ", $this_genome_db->name, " has ", scalar(@$all_genes), " $biotype genes\n";
    foreach my $this_gene (@$all_genes) {

      my $this_species_name = $this_genome_db->name;
      my $this_gene_slice = $this_gene->slice->sub_Slice($this_gene->start, $this_gene->end);
      my $num_of_synteny_regions = 0;
      my $max_overlap = 0;
      my $total_overlap = 0;
      my $max_overlap_dnafrag_region = undef;
      my $last_end = undef;
      foreach my $this_dnafrag_region (@{$all_dnafrag_regions->{$this_species_name}->{$this_gene_slice->seq_region_name}}) {
        my $overlap_slice;
        next if ($this_dnafrag_region->dnafrag_end < $this_gene_slice->start);
        last if ($this_dnafrag_region->dnafrag_start > $this_gene_slice->end);
        $overlap_slice = get_overlap($this_gene_slice, $this_dnafrag_region->slice)
            if ($this_dnafrag_region->dnafrag->genome_db->name eq $this_genome_db->name);
        my $overlap_length;
        if ($overlap_slice) {
          $num_of_synteny_regions++;
          if ($last_end) {
            if ($overlap_slice->end <= $last_end) {
              next;
            }
            if ($overlap_slice->start <= $last_end) {
              # Test for the max_overlap before truncating the overlap slice
              $overlap_length = $overlap_slice->length;
              if ($overlap_length > $max_overlap) {
                $max_overlap_dnafrag_region = $this_dnafrag_region;
                $max_overlap = $overlap_length;
              }
              # Truncate the overlap slice, so no bp will be counted twice in the total
              $overlap_slice = $overlap_slice->sub_Slice($last_end - $overlap_slice->start + 2, $overlap_slice->length);
            }
          }
          $overlap_length = $overlap_slice->length;
          if ($overlap_length > $max_overlap) {
            $max_overlap_dnafrag_region = $this_dnafrag_region;
            $max_overlap = $overlap_length;
          }
          $total_overlap += $overlap_length;
          $last_end = $overlap_slice->end;
        }
      }
      if ($num_of_synteny_regions == 0) {
        $uncovered++;
      } elsif ($max_overlap == $this_gene_slice->length) {
        $fully++;
      } elsif ($total_overlap == $this_gene_slice->length) {
        $broken++;
      } else {
        $partially++;
      }
    }
    if (scalar(@$all_genes)) {
      printf "  $fully (%.1f%%) fully covered; $broken (%.1f%%) broken; $partially (%.1f%%) partially; $uncovered (%.1f%%) missing\n",
          ($fully*100/scalar(@$all_genes)), ($broken*100/scalar(@$all_genes)),
          ($partially*100/scalar(@$all_genes)), ($uncovered*100/scalar(@$all_genes));
    } else {
      printf "  $fully (%.1f%%) fully covered; $broken (%.1f%%) broken; $partially (%.1f%%) partially; $uncovered (%.1f%%) missing\n",
          0, 0, 0, 0;
    }
  #   for (my $i = 0; $i < @$totals; $i++) {
  #     print " $i\t", ($totals->[$i] or 0), "\t", ($maxs->[$i] or 0), "\t", ($sizes->[$i] or 0), "\n";
  #   }
  #   print join("\n", map {$_+=0} @$maxs), "\n";
  }
}


sub get_overlap {
  my ($slice1, $slice2) = @_;

  return undef if ($slice1->coord_system_name ne $slice2->coord_system_name);
  return undef if ($slice1->seq_region_name ne $slice2->seq_region_name);
  if ($slice1->start < $slice2->end and $slice2->start < $slice1->end) {
    my $start = ($slice1->start < $slice2->start)?$slice2->start:$slice1->start;
    my $end = ($slice1->end < $slice2->end)?$slice1->end:$slice2->end;
    return $slice1->seq_region_Slice->sub_Slice($start, $end);
  }
  return undef;
}


sub get_all_genes {
  my ($table_name, $species) = @_;

  my $sth;
  if ($species) {
    $sth = $compara_dba->dbc->prepare(qq!SELECT * FROM $table_name WHERE species_name = ?!);
    $sth->execute($species);
  } else {
    $sth = $compara_dba->dbc->prepare(qq!SELECT * FROM $table_name!);
    $sth->execute();
  }
  return $sth->fetchall_arrayref({});
}


=head2 get_MethodLinkSpeciesSet

  Arg[1]      : Bio::EnsEMBL::Compara::DBSQL::DBAdaptor $compara_dba
  Arg[1]      : int $method_link_species_set_id
  Description : Fetches the Bio::EnsEMBL::Compara::MethodLinkSpeciesSet
                object corresponding to this $method_link_species_set_id.
  Returns     : Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object
  Exceptions  :

=cut

sub get_MethodLinkSpeciesSet {
  my ($compara_dba, $method_link_species_set_id) = @_;

  # Getting Bio::EnsEMBL::Compara::MethodLinkSpeciesSet obejct
  my $method_link_species_set_adaptor = $compara_dba->get_MethodLinkSpeciesSetAdaptor;

  my $method_link_species_set =
      $method_link_species_set_adaptor->fetch_by_dbID($method_link_species_set_id);
# print "Using ", $method_link_species_set->name, "\n";

  throw("The database do not contain MethodLinkSpeciesSet with dbID $method_link_species_set_id!")
      if (!$method_link_species_set);

  return $method_link_species_set;
}


=head2 get_all_DnaFragRegions

  Arg[1]      : Bio::EnsEMBL::Compara::DBSQL::DBAdaptor $compara_dba
  Arg[1]      : Bio::EnsEMBL::Compara::MethodLinkSpeciesSet $method_link_species_set
  Description : Fetches all the Bio::EnsEMBL::Compara::SyntenyRegion
                objects for this $method_link_species_set in this DB.
  Returns     : listref of Bio::EnsEMBL::Compara::SyntenyRegion objects
  Exceptions  :

=cut

sub get_all_DnaFragRegions {
  my ($all_synteny_region) = @_;
  my $all_dnafrag_regions;

  foreach my $this_synteny_region (@$all_synteny_regions) {
    foreach my $this_dnafrag_region (@{$this_synteny_region->get_all_DnaFragRegions()}) {
      my $dnafrag = $this_dnafrag_region->dnafrag;
      my $species_name = $dnafrag->genome_db->name;
      my $dnafrag_name = $dnafrag->name;
      push(@{$all_dnafrag_regions->{$species_name}->{$dnafrag_name}}, $this_dnafrag_region);
    }
  }
  foreach my $this_species (keys %$all_dnafrag_regions) {
    foreach my $this_chr_name (keys %{$all_dnafrag_regions->{$this_species}}) {
      @{$all_dnafrag_regions->{$this_species}->{$this_chr_name}} = sort {$a->dnafrag_start <=> $b->dnafrag_start}
          @{$all_dnafrag_regions->{$this_species}->{$this_chr_name}};
    }
  }

  return $all_dnafrag_regions;
}


=head2 fetch_all_SyntenyRegions

  Arg[1]      : Bio::EnsEMBL::Compara::DBSQL::DBAdaptor $compara_dba
  Arg[1]      : Bio::EnsEMBL::Compara::MethodLinkSpeciesSet $method_link_species_set
  Description : Fetches all the Bio::EnsEMBL::Compara::SyntenyRegion
                objects for this $method_link_species_set in this DB.
  Returns     : listref of Bio::EnsEMBL::Compara::SyntenyRegion objects
  Exceptions  :

=cut

sub fetch_all_SyntenyRegions {
  my ($compara_dba, $method_link_species_set) = @_;

  # Getting Bio::EnsEMBL::Compara::MethodLinkSpeciesSet obejct
  my $method_link_species_set_adaptor = $compara_dba->get_MethodLinkSpeciesSetAdaptor;

  my $method_link_species_set =
      $method_link_species_set_adaptor->fetch_by_dbID($method_link_species_set_id);
# print "Using ", $method_link_species_set->name, "\n";

  throw("The database do not contain MethodLinkSpeciesSet with dbID $method_link_species_set_id!")
      if (!$method_link_species_set);

  my $synteny_region_adaptor = $compara_dba->get_SyntenyRegionAdaptor;
  my $all_synteny_regions = $synteny_region_adaptor->fetch_all_by_MethodLinkSpeciesSet(
      $method_link_species_set);

  return $all_synteny_regions;
}


=head2 print_genomic_coverage

  Arg[1]      : listref Bio::EnsEMBL::Compara::SyntenyRegion $all_synteny_regions
  Arg[1]      : [optional] int $verbose
  Description :
  Returns     :
  Exceptions  :

=cut

sub print_genomic_coverage {
  my ($all_dnafrag_regions, $quick) = @_;

  ## Get the length for all the DnaFrags (well, only if they have been used in a SyntenyRegion)
  my $dnafrag_length;
  foreach my $species_name (keys %{$all_dnafrag_regions}) {
    foreach my $this_dnafrag_name (keys %{$all_dnafrag_regions->{$species_name}}) {
      my $dnafrag = $all_dnafrag_regions->{$species_name}->{$this_dnafrag_name}->[0]->dnafrag;
      $dnafrag_length->{$species_name}->{$this_dnafrag_name} = $dnafrag->length;
    }
  }

  ## Loop through all the species
  while (my ($species_name, $this_ref) = each %$all_dnafrag_regions) {
    if ($quick) {
      print "Species name              : CHR --  \%COV -- \%UNCV\n";
    } else {
      print "Species name              : CHR --  \%COV --  \%GAP --  \%REP -- \%REST\n";
    }
    my $coverage_per_species = 0;
    my $ass_gap_non_covered_per_species = 0;
    my $repeat_non_covered_per_species = 0;
    my $total_length = 0;
    my $longest_uncovered_regions;

    ## Loop through all the DnaFrags for this species
    foreach my $dnafrag_name (sort {($a=~/^\d+$/ and $b=~/^\d+$/ and $a<=>$b) or ($a cmp $b)} keys %$this_ref) {
      my $all_regions = $this_ref->{$dnafrag_name};
#       @$all_regions = sort {$a->dnafrag_start <=> $b->dnafrag_start} @$all_regions;
      my $coverage = 0;
      my $ass_gap_non_covered = 0;
      my $repeat_non_covered = 0;
      my $last_end = 0;
      my $this_slice = $all_regions->[0]->dnafrag->slice;
      my $coord_system_name = $all_regions->[0]->dnafrag->coord_system_name;
      foreach my $this_region (@$all_regions) {
        if ($this_region->dnafrag_start <= $last_end) {
          ## Overlaps previous region
          if ($this_region->dnafrag_end > $last_end) {
            $coverage += $this_region->dnafrag_end - $last_end; # $this_region->dnafrag_end - ($last_end + 1) + 1
            $last_end = $this_region->dnafrag_end;
          }
        } else {
          if (!$quick and $coord_system_name eq "chromosome") {
            if ($last_end + 1 < $this_region->dnafrag_start - 1) {
              my $slice = $this_slice->sub_Slice($last_end + 1, $this_region->dnafrag_start - 1);
              $ass_gap_non_covered += get_length_of_assembly_gaps($slice);
              $repeat_non_covered += get_length_of_repeats($slice);
            }
          }

          $coverage += $this_region->dnafrag_end - $this_region->dnafrag_start + 1;
          $last_end = $this_region->dnafrag_end;
        }
      }
      # Get assembly gap and repeats on the last bit of the chromosome
      if (!$quick and $coord_system_name eq "chromosome") {
        if ($last_end + 1 < $dnafrag_length->{$species_name}->{$dnafrag_name}) {
          my $slice = $this_slice->sub_Slice($last_end + 1, $dnafrag_length->{$species_name}->{$dnafrag_name});
          $ass_gap_non_covered += get_length_of_assembly_gaps($slice);
          $repeat_non_covered += get_length_of_repeats($slice);
        }
      }

      ## Print values for this DnaFrag
      if ($quick) {
        printf "%-25s : %3s -- %5.2f -- %5.2f\n",
            $species_name,
            $dnafrag_name,
            ($coverage*100/$dnafrag_length->{$species_name}->{$dnafrag_name}),
            (($dnafrag_length->{$species_name}->{$dnafrag_name} - $coverage)*100/
                $dnafrag_length->{$species_name}->{$dnafrag_name});
      } else {
        printf "%-25s : %3s -- %5.2f -- %5.2f -- %5.2f -- %5.2f\n",
            $species_name,
            $dnafrag_name,
            ($coverage*100/$dnafrag_length->{$species_name}->{$dnafrag_name}),
            ($ass_gap_non_covered*100/$dnafrag_length->{$species_name}->{$dnafrag_name}),
            ($repeat_non_covered*100/$dnafrag_length->{$species_name}->{$dnafrag_name}),
            (($dnafrag_length->{$species_name}->{$dnafrag_name} - $coverage - $ass_gap_non_covered -
                $repeat_non_covered)*100/$dnafrag_length->{$species_name}->{$dnafrag_name});
      }
      $coverage_per_species += $coverage;
      $ass_gap_non_covered_per_species += $ass_gap_non_covered;
      $repeat_non_covered_per_species += $repeat_non_covered;
      $total_length += $dnafrag_length->{$species_name}->{$dnafrag_name};
    }

    ## Print final line, summary for this species
    if ($quick) {
      printf "%-25s : ALL -- %5.2f -- %5.2f\n\n",
          $species_name,
          ($coverage_per_species*100/$total_length),
          (($total_length - $coverage_per_species - $ass_gap_non_covered_per_species -
              $repeat_non_covered_per_species)*100/$total_length);
    } else {
      printf "%-25s : ALL -- %5.2f -- %5.2f -- %5.2f -- %5.2f\n\n",
          $species_name,
          ($coverage_per_species*100/$total_length),
          ($ass_gap_non_covered_per_species*100/$total_length),
          ($repeat_non_covered_per_species*100/$total_length),
          (($total_length - $coverage_per_species - $ass_gap_non_covered_per_species -
              $repeat_non_covered_per_species)*100/$total_length);
    }

  }
}


=head2 print_N50_stats

  Arg[1]      : listref Bio::EnsEMBL::Compara::SyntenyRegion $all_synteny_regions
  Arg[1]      : [optional] int $verbose
  Description :
  Returns     :
  Exceptions  :

=cut

sub print_N50_stats {
  my ($all_synteny_regions) = @_;

  while (my ($species_name, $this_ref) = each %$all_dnafrag_regions) {
    my $lengths;
    my $sum_length = 0;

    ## Loop through all the DnaFrags for this species
    foreach my $all_regions (values %$this_ref) {
      my $coverage = 0;
      foreach my $this_region (@$all_regions) {
        my $this_length = $this_region->dnafrag_end - $this_region->dnafrag_start + 1; 
        push(@$lengths, $this_length);
        $sum_length += $this_length;
      }
    }

    my $partial_sum = 0;
    my $n50;
    foreach my $this_length (sort {$b <=> $a} @$lengths) {
      $partial_sum += $this_length;
      if ($partial_sum >= $sum_length / 2) {
        $n50 = $this_length;
        last;
      }
    }
    ## Print result line
    printf "%-25s : N50 = %10d\n", $species_name, $n50;
  }
  print "\n";
}


=head2 print_duplication_stats

  Arg[1]      : listref Bio::EnsEMBL::Compara::SyntenyRegion $all_synteny_regions
  Arg[1]      : [optional] int $verbose
  Description :
  Returns     :
  Exceptions  :

=cut

sub print_duplication_stats {
  my ($all_synteny_regions, $species_tree, $verbose) = @_;

  my $duplications = [];
  foreach my $this_synteny_region (@$all_synteny_regions) {
    my $species;
    foreach my $this_dnafrag_region (@{$this_synteny_region->get_all_DnaFragRegions()}) {
      my $species_name = $this_dnafrag_region->genome_db->name;
      $species->{$species_name}++;
    }
    foreach my $cardinality (values %$species) {
      if ($cardinality > 1) {
        push(@$duplications, [$this_synteny_region, $species]);
        last;
      }
    }
  }
  print "Found ", scalar(@$duplications), " duplications\n";

  my $duplications_per_species;
  my $ancient_duplications = [];
  my $recent_count = 0;
  my $single_count = 0;
  foreach my $object (@$duplications) {
    my ($this_synteny_region, $species) = @$object;
    my $duplicated_genome_dbs = {};

    foreach my $this_dnafrag_region (@{$this_synteny_region->get_all_DnaFragRegions()}) {
      my $species_name = $this_dnafrag_region->genome_db->name;
      if ($species->{$species_name} > 1) {
        $duplicated_genome_dbs->{$species_name} = $this_dnafrag_region->genome_db;
        $duplications_per_species->{$species_name} += $this_dnafrag_region->length;
      }
    }

    if (scalar(values(%$species)) > 1) {
      my @dup_genomes = values %$duplicated_genome_dbs;
      if (@dup_genomes > 1) {
        my $ancestor = find_shared_ancestor($dup_genomes[0]->taxon, $dup_genomes[1]->taxon);
        foreach my $this_genome (@dup_genomes[2..$#dup_genomes]) {
          $ancestor = find_shared_ancestor($ancestor, $this_genome->taxon);
        }
        if ($ancestor->node_id == $species_tree->node_id) {
#           print "DUP GENOME_DBS: ", join(" -- ", sort map {$_->name} @dup_genomes), "\n";
          push(@$ancient_duplications, $this_synteny_region);
        }
      } else {
        $recent_count++;
      }
    } else {
      $single_count++;
    }
  }

  foreach my $this_species (sort keys %$duplications_per_species) {
    print " Duplications on $this_species span ", $duplications_per_species->{$this_species}, " bp\n";
  }
  print "\n";

  print "Found ", scalar(@$ancient_duplications), " ancient duplications\n";
  print "Found $recent_count recent duplications\n";
  print "Found $single_count species-specific duplicated regions\n";

  foreach my $this_ancient_duplication (@$ancient_duplications) {
    foreach my $this_dnafrag_region (sort {$a->dnafrag->genome_db->name cmp $b->dnafrag->genome_db->name ||
        $a->dnafrag->name cmp $b->dnafrag->name || $a->dnafrag_start <=> $b->dnafrag_start}
        @{$this_ancient_duplication->get_all_DnaFragRegions()}) {
      my $species_name = $this_dnafrag_region->genome_db->name;
      $species_name =~ s/ /_/g;
      print "http://www.ensembl.org/$species_name/contigview?region=", $this_dnafrag_region->dnafrag->name,
          ";start=", $this_dnafrag_region->dnafrag_start, ";end=", $this_dnafrag_region->dnafrag_end,
          " [", (($this_dnafrag_region->dnafrag_strand == -1)?"-":"+"), "](", $this_dnafrag_region->length, ")\n";
    }
    print "\n";
  }
}


=head2 get_length_of_assembly_gaps

  Arg[1]      : Bio::EnsEMBL::Slice $slice
  Description : This methods looks for gaps in the assembly lying on this
                particular slice. The total length of all the assembly
                gaps is returned as an integer.
  Returns     : integer
  Exceptions  :

=cut

sub get_length_of_assembly_gaps {
  my ($slice) = @_;
  my $assembly_gap_length = 0;

  my $proj_segments = $slice->project("seqlevel");
  my $last_end = 0;
  foreach my $this_proj_segment (@$proj_segments) {
    ## last_end refers to the end of previous segment. ($last_end + 1) would be
    ## the next position. Threfore the difference between actual from_start and
    ## ($last_end + 1) is the length of the gap (0 if no gap)
    $assembly_gap_length += $this_proj_segment->from_start - ($last_end + 1);
    $last_end = $this_proj_segment->from_end;
  }

  return $assembly_gap_length;
}


=head2 get_length_of_repeats

  Arg[1]      : Bio::EnsEMBL::Slice $slice
  Description : This methods looks for repeat features lying on this
                particular slice. The total length of all these repeats
                is returned as an integer. Note that this method is aware
                of overlapping repeats, ie overlapping repeats will be
                counted only once.
  Returns     : integer
  Exceptions  :

=cut

sub get_length_of_repeats {
  my ($slice) = @_;
  my $repeat_length = 0;

  my $repeat_features = $slice->get_all_RepeatFeatures;
  my $last_end = 0;
  foreach my $this_repeat_feature (sort {$a->start <=> $b->start} @$repeat_features) {
    if ($this_repeat_feature->start <= $last_end) {
      ## Overlaps previous feature
      if ($this_repeat_feature->end > $last_end) {
        $repeat_length += $this_repeat_feature->end - $last_end; # $this_repeat_feature->dnafrag_end - ($last_end + 1) + 1
        $last_end = $this_repeat_feature->end;
      }
    } else {
      $repeat_length += $this_repeat_feature->end - $this_repeat_feature->start + 1;
      $last_end = $this_repeat_feature->end;
    }
  }

  return $repeat_length;
}


sub get_tree_by_GenomeDBs {
  my ($genome_dbs) = @_;

  my $first_genome_db = shift @$genome_dbs;
  my $node = $first_genome_db->taxon;
  $node->no_autoload_children;
  my $root = $node->root;

  foreach my $genome_db (@$genome_dbs) {
    my $node = $genome_db->taxon;
    unless (defined $node) {
      print STDERR $genome_db->name, " not in the database\n";
      next;
    }
    $node->no_autoload_children;
    $root->merge_node_via_shared_ancestor($node);
  }

  $root = $root->minimize_tree;

  return $root;
}

sub find_shared_ancestor {
  my $first_node = shift;
  my $second_node = shift;

  my $exit1 = 0;
  while (!$exit1) {
    my $tmp_node = $second_node;
    my $exit2 = 0;
    while (!$exit2) {
      if ($tmp_node->node_id == $first_node->node_id) {
        return $first_node;
      } elsif ($tmp_node->parent) {
        $tmp_node = $tmp_node->parent;
      } else {
        $exit2 = 1;
      }
    }
    if ($first_node->parent) {
      $first_node = $first_node->parent;
    } else {
      $exit1;
    }
  }

  return undef;
}

=head1 OLD RESULTS

=head2 Case A

# max-path-dissimilarity: 0 # max-ratio: off # simplify-graph: 0

perl ../../tools/assess_graph.pl --alignment_type ENREDO-A --set_of_species human:mouse

+ Mus musculus has 24442 protein_coding genes
B<10765 (44.0%) fully covered; 7261 (29.7%) partially/broken; 6416 (26.2%) missing>

+ Homo sapiens has 22810 protein_coding genes
B<10126 (44.4%) fully covered; 7640 (33.5%) partially/broken; 5044 (22.1%) missing>

=head2 Case B

# max-path-dissimilarity: 3 # max-ratio: off # simplify-graph: 0

perl ../../tools/assess_graph.pl --alignment_type ENREDO-B --set_of_species human:mouse

+ Mus musculus has 24442 protein_coding genes
B<14542 (59.5%) fully covered; 4547 (18.6%) partially/broken; 5353 (21.9%) missing>

+ Homo sapiens has 22810 protein_coding genes
B<13776 (60.4%) fully covered; 4838 (21.2%) partially/broken; 4196 (18.4%) missing>

=head2 Case C

# max-path-dissimilarity: 3 # max-ratio: off # simplify-graph: 1

$ perl ../../tools/assess_graph.pl --alignment_type ENREDO-C --set_of_species human:mouse

+ Mus musculus has 24442 protein_coding genes
B<15895 (65.0%) fully covered; 3294 (13.5%) partially/broken; 5253 (21.5%) missing>

+ Homo sapiens has 22810 protein_coding genes
B<15227 (66.8%) fully covered; 3565 (15.6%) partially/broken; 4018 (17.6%) missing>

=head2 Case D

# max-path-dissimilarity: 3 # max-ratio: off # simplify-graph: 2

$ perl ../../tools/assess_graph.pl --alignment_type ENREDO-D --set_of_species human:mouse

+ Mus musculus has 24442 protein_coding genes
B<15932 (65.2%) fully covered; 3272 (13.4%) partially/broken; 5238 (21.4%) missing>

+ Homo sapiens has 22810 protein_coding genes
B<15269 (66.9%) fully covered; 3558 (15.6%) partially/broken; 3983 (17.5%) missing>

=head2 Case E

# max-path-dissimilarity: 3 # max-ratio: off # simplify-graph: 3

+ Mus musculus has 24442 protein_coding genes
B<16642 (68.1%) fully covered; 2753 (11.3%) partially/broken; 5047 (20.6%) missing>

+ Homo sapiens has 22810 protein_coding genes
B<15967 (70.0%) fully covered; 3000 (13.2%) partially/broken; 3843 (16.8%) missing>

=head2 Case F

# max-path-dissimilarity: 3 # max-ratio: off # simplify-graph: 4

+ Mus musculus has 24442 protein_coding genes
B<16643 (68.1%) fully covered; 2752 (11.3%) partially/broken; 5047 (20.6%) missing>

+ Homo sapiens has 22810 protein_coding genes
B<15968 (70.0%) fully covered; 3000 (13.2%) partially/broken; 3842 (16.8%) missing>

=head2 Case G

# max-path-dissimilarity: 3 # max-ratio: 3.0 # simplify-graph: 4

+ Mus musculus has 24442 protein_coding genes
B<16773 (68.6%) fully covered; 2666 (10.9%) partially/broken; 5003 (20.5%) missing>

+ Homo sapiens has 22810 protein_coding genes
B<16092 (70.5%) fully covered; 2916 (12.8%) partially/broken; 3802 (16.7%) missing>

=head2 Case H

# max-path-dissimilarity: 3 # max-ratio: 2.0 # simplify-graph: 4

+ Mus musculus has 24442 protein_coding genes
B<16792 (68.7%) fully covered; 2642 (10.8%) partially/broken; 5008 (20.5%) missing>

+ Homo sapiens has 22810 protein_coding genes
B<16107 (70.6%) fully covered; 2895 (12.7%) partially/broken; 3808 (16.7%) missing>

=head2 Case I

# max-path-dissimilarity: 5 # max-ratio: 2.0 # simplify-graph: 1

+ Mus musculus has 24442 protein_coding genes
B<16091 (65.8%) fully covered; 3184 (13.0%) partially/broken; 5167 (21.1%) missing>

+ Homo sapiens has 22810 protein_coding genes
B<15410 (67.6%) fully covered; 3452 (15.1%) partially/broken; 3948 (17.3%) missing>

=head2 Case J

# max-path-dissimilarity: 5 # max-ratio: 2.0 # simplify-graph: 4

+ Mus musculus has 24442 protein_coding genes
B<17217 (70.4%) fully covered; 2390 (9.8%) partially/broken; 4835 (19.8%) missing>

+ Homo sapiens has 22810 protein_coding genes
B<16474 (72.2%) fully covered; 2641 (11.6%) partially/broken; 3695 (16.2%) missing>

=cut
