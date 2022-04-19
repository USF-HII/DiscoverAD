# DiscoverAD

Discover Atypical Diabetes (DiscoverAD) was developed for identifying and clustering endotypes of atypical diabetes. DiscoverAD is a data mining framework with a two-step filtering process to first exclude participants who meet definitions of typical type 1 diabetes (T1D) or type 2 diabetes (T2D) and then include participants with certain pre-specified AD characteristics. This is followed by robust and unsupervised cluster analysis to discover novel endotypes of AD within the filtered group. We purposefully developed DiscoverAD to permit flexibility and efficiency so it can be applicable for various clinical settings with different types of large cohort datasets.

DiscoverAD is written in the R statistical programming language. A tutorial is available [here](DiscoverAD_tutorial.md) in markdown format or as plain [R Code](DiscoverAD_tutorial.R). Additionally, we've included sample files for the tutorial:
1.  [Example dataset format](DiscoverAD_tutorial_files/sample_input_files/DiscoverAD_sample_data.tsv) - an example dataset to demonstrate the format used for DiscoverAD.
2.  [Example variable type defintion file](DiscoverAD_tutorial_files/sample_input_files/variable_type.tsv) - a file including the statistical [data type](https://github.com/USF-HII/DiscoverAD/blob/master/DiscoverAD_tutorial.md#input-data) of each variable included in the example dataset.
