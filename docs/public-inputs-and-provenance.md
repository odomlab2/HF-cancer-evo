# Public inputs and provenance

The advertised public command reads these tracked input classes:

- manuscript plot-data CSV files beneath `Figures/`;
- the frozen human-gene incidence table and permutation draws beneath
  `Figures/Figure 4/`;
- compact submitted-statistics tables beneath `publication_inputs/`.

The compact statistics tables retain publication-critical values that were
already tracked in the submitted-analysis evidence. Private source-object paths
were removed. They are validation inputs, not replacements for undistributed
primary sequencing data.

The human incidence and clinical plot-data tables are frozen cBioPortal-derived
submission inputs. The active analysis names the four public cSCC studies and
the public cBioPortal API endpoint. MSigDB-named files are frozen annotation
exports used for manuscript plot data, not redistributed gene-set collections.
No live service is contacted by the advertised command.

Primary sequencing files, raw QC exports, reference FASTA files, VEP caches,
generated RDS objects, and reference CDS objects are external inputs. They are
not distributed in this snapshot.
