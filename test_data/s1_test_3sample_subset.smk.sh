snakemake -s Snakefile4test --configfile test_data/config4test_sc.3sample_subset.yaml --forceall -j 1 test_sc_rna_default
snakemake -s Snakefile4test --configfile test_data/config4test_sc.3sample_subset.yaml --forceall -j 1 test_sc_rna_exclude_introns 
snakemake -s Snakefile4test --configfile test_data/config4test_sc.3sample_subset.yaml --forceall -j 1 test_sc_rna_expect_cell
snakemake -s Snakefile4test --configfile test_data/config4test_sc.3sample_subset.yaml --forceall -j 1 test_sc_rna_force_cell
snakemake -s Snakefile4test --configfile test_data/config4test_sc.3sample_subset.yaml --forceall -j 1 test_sc_multi_default 
snakemake -s Snakefile4test --configfile test_data/config4test_sc.3sample_subset.yaml --forceall -j 1  test_sc_multi_expect_cell 
snakemake -s Snakefile4test --configfile test_data/config4test_sc.3sample_subset.yaml --forceall -j 1  test_sc_multi_force_cell
snakemake -s Snakefile4test --configfile test_data/config4test_sc.3sample_subset.yaml --forceall -j 1 test_sc_vdj_default
snakemake -s Snakefile4test --configfile test_data/config4test_sc.3sample_subset.yaml --forceall -j 1 test_sc_vdj_chain
snakemake -s Snakefile4test --configfile test_data/config4test_sc.3sample_subset.yaml --forceall -j 1 test_sc_multiome_default 
snakemake -s Snakefile4test --configfile test_data/config4test_sc.3sample_subset.yaml --forceall -j 1 test_sc_multiome_exclude_exons
snakemake -s Snakefile4test --configfile test_data/config4test_sc.3sample_subset.yaml --forceall -j 1 test_sc_spatial_default

