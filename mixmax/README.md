#mixmax

* `create_demultiplexing_scheme.py` takes the number of samples and the maximal pool size as input and outputs a demultiplexin scheme.

  * Usage: `python create_demultiplexing_scheme.py -- robust {True,False} --maximal_pool_size {k} --n_samples {n} --input_dir {dir} --output {file}`

  * `--robust` is a boolean value that specifies if the output multiplexing scheme should be robust.
  * `--maximal_pool_size` specifies that no pool should contain more than `k`samples.
  * `--n_samples` defines the size of the cohort. 
  * `--output` specifies that the multiplexing scheme should be output to `file`. 

* `match_samples.py`: takes the genotypes of demultiplexed datasets and outputs a file containing the sample identification information.
  * Usage: `python --tsv_files {file1.tsv file2.tsv} --pool_scheme {scheme.yaml} --output_plot {plot.png} --output {sample_assignment.yaml}  --robust {True,False}
  
  * `--tsv-files` contains all demultplexed genotypes, one tsv-file per pool. It expects rows to be datasets and columns to be genomic position. The entries are values between 0 and 1 (the expected variant allele frequency in the dataset).
  * `pool_scheme` expect a yaml-file as output by `create_demultiplexing_scheme.py`.
  * `--output_plot` takes the path to which a diagnostic heatmap will be saves to.
  * `--output` expect the path to the `sample_assignment.yaml` file which specifies the sample identification in human readable format.
  * `--robust` must be specified as in `create_demultiplexing_scheme.py`.