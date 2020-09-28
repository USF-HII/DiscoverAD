# Strip the R sections out of the tutorial and write as a simple R script for convenience

all:
	sed -n '/^``` r$$/,/^```$$/p' DiscoverAD_tutorial.md | sed 's/^```.*//' > DiscoverAD_tutorial.R

