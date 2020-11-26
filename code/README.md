### Table of Contents  
[Main settings](#main-settings)  
[Additional settings](#additional-settings)  
[Input data](#input-data)  
[Launch the pipeline](#launch-the-pipeline)  
[Pipeline content](#pipeline-content)  
[Outputs](#outputs)  
[Dependencies](#depedencies)  
[Warning](#warning)  
[Example](#example)  
[Script dependencies](#script-dependencies)  

### Main settings
([Back to top](#table-of-contents))

Set main settings in <em>main_settings.R</em> (cf. template file: <em>TEMPLATE_main_settings.R</em>).

Do not change the name of this file.

The current <em>main_settings.R</em> file was used to run the example as explained below.


### Additional settings
([Back to top](#table-of-contents))

Set additional (e.g. dataset-specific) settings in another file (cf. template file: <em>TEMPLATE_setting_file.R</em>)

The name of this file can be user-defined.

Settings from this file will overwrite settings from <em>main_settings.R</em> (hence can be used for this purpose).


### Input data
([Back to top](#table-of-contents))

The expected format for the input data is described <a href="https://raw.githack.com/marzuf/MANUSCRIPT_FIGURES/master/code/data/input_data_desc.html">here</a>. 

<!-- The expected format of the input data is described <a href="https://htmlpreview.github.io/?https://github.com/marzuf/MANUSCRIPT_FIGURES/blob/master/code/data/input_data_desc.html">here</a>. -->
<!-- TOC not rendering with html preview -->


### Launch the pipeline
([Back to top](#table-of-contents))

```{bash}
 ./run_pipeline.sh <setting_file> <stepNbr> <stepNbr> <...>
```

`setting_file` should be path to the additional settings file.

`stepNbr` indicates the step(s) of the pipeline that should be run (should be what comes before the "_" in the script name, e.g. `1` or `5fc`)

### Pipeline content
([Back to top](#table-of-contents))

The pipeline is designed to run the steps one after the other (some [depedencies](#script-dependencies) exist among the scripts, e.g. the 11th requires the 10th and 9th to have been run).

The names of the scripts explicitly describe to which step they correspond, namely:
<ul>
<li><em>1_prepGeneData.R</em>: prepare some data for the rest of the pipeline</li>
<li><em>2_runGeneDE.R</em>: gene-level differential expression analysis (needed for log2-fold-change (logFC) values)</li>
<li><em>3_runMeanTADLogFC.R</em>: prepare domain-level average logFC</li>
<li><em>4_runMeanTADCorr.R</em>: prepare domain-level average correlation</li>
<li><em>5fc_runPermutationsMedian.R</em>: run permutation for logFC values (! can take some time - especially for saving the file - (and disk space) depending on the number of permutations !; from our experience (40 CPU): ~20 min for 1000 permutations, ~2 h for 100'000)</li>
<li><em>5corr_runPermutationsCorr.R</em>: run permutation for correlation values</li>
<li><em>6_runPermutationsMeanLogFC.R</em>: compute domain-level average FC for the permutation data (! take ca. ~2 min for 1000 permutations, ~2 h for 100'000)</li>
<li><em>7_runPermutationsMeanTADCorr.R</em>:  compute domain-level average correlation for the permutation data</li>
<li><em>8_runEmpPvalMeanTADLogFC.R</em>: compute empirical p-values for the logFC</li>
<li><em>9_runEmpPvalMeanTADCorr.R</em>: compute empirical p-values for correlation</li>
<li><em>10_runEmpPvalCombined.R</em>: combine empirical p-values for FC and correlation</li>
<li><em>11_runFCC.R</em>: fold-change concordance (FCC) for the observed data</li>
<li><em>12_runPermutationsFCC.R</em>: FCC for the permutation data (data from from 5fc; ! take ca. ~2 min for 1000 permutations, ~2h15 for 100'000)</li>
<li><em>13_runFCCcumSumAUC.R</em>: obs/permut. AUC ratio of FCC cumsum curves</li>
</ul>

### Outputs
([Back to top](#table-of-contents))

In the folder indicated in the settings, each script mainly outputs RData files in a folder with the same name as the script (e.g. files written in <em>1_prepGeneData</em> for <em>1_prepGeneData.R</em>).


### Dependencies
([Back to top](#table-of-contents))


[pigz](https://zlib.net/pigz) is required for fast saving in RData files in the scripts 5fc, 6 and 12 (we used version 2.4).


R packages used in the pipeline: limma and edge  (for gene-level differential expression analysis in step 2), foreach, doMC, dplyr, tools, flux.


The pipeline has been tested in July 2020 with R version 3.6.0.


### **WARNING** 
([Back to top](#table-of-contents))

In step 9 (also impacting step 10), there are two possibilities for retrieving correlation permutation data to be set in the setting files:

<ul>
<li>if <code>all_permutCorr_data</code> is the path to a folder:</li>
<ul>
<li>the correlation values from permutation data are loaded from files in <code>all_permutCorr_data</code> that (recursively) match the pattern <code>corrMatchPattern</code> (currently set to <code>"meanCorr_sample_around_TADs_sameNbr.RData"</code> (output of step 7) in the main settings but can be overwritten in the additional settings)</li>
<li>if <code>refineMatchPattern</code> is provided, this pattern is used for a second pattern matching to refine file retrieval</li>
<li>if <code>corrDiscardPattern</code> is provided in setting files, files that match that pattern are discarded</li>
<li>if <code>nbrCorrPermutCheck</code> is provided, it is checked that as many as <code>nbrCorrPermutCheck</code> files have been loaded</li>
<li>it is expected that retrieved files correspond to the output of script 7: each file should be loadable (RData) and contain a list of lists storing a correlation value in <code>my_list[[idx]][["meanCorr"]]</code> with as many <code>idx</code> as domains
<li>ouput files from steps 9 and 10 will be prefixed with <em>fromFolder_</em>
</ul>
<li>if <code>all_permutCorr_data</code> is the path to a file:</li>
<ul>
<li><code>all_permutCorr_file</code> should provide the path to correlation values for the permutation data (the file <em>data/all_sample_corrValues.RData</em> provides the values from our permuation data); the data should be a RData file containing a list/vector of permutation correlation values</li>
<li>ouput files from steps 9 and 10 will be prefixed with <em>fromFile_</em>
</ul>
</ul>


### Example
([Back to top](#table-of-contents))

Example for running the full pipeline with our LUAD data (using provided data for intra-domain correlation, cf. setting files):

```{bash}
 ./example_run_pipeline.sh example_ENCSR489OCU_NCI-H460_40kb_run_settings_TCGAluad_norm_luad.R 1 2 3 4 5fc 5corr 6 7 8 9 10 11 12 13
```

(ca. 55 minutes to run the full the pipeline with 40 CPU and 1000 permutations, 8h for 100'000 permutations)


### Script dependencies
([Back to top](#table-of-contents))

| Script        | Dependencies      |
| ------------- |-------------------|
| 1             | -                 |
| 2             | 1                 |
| 3             | 1, 2              |
| 4             | 1                 |
| 5fc           | 1                 |
| 5corr         | 1                 |
| 6             | 1, 2 , 5fc        |
| 7             | 1, 5corr          |
| 8             | 3, 6              |
| 9             | 4, 7              |
| 10            | 8, 9              |
| 11            | 1, 2              |
| 12            | 1, 2, 5fc         |
| 13            | 11, 12            |



