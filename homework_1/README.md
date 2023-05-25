# Differential gene expression #

This script `diffexpress.py` has several functions that calculate the difference between gene expressions of two cell lines.

## Main script function ##

**The main function of the script** is `diffexpr` that takes 3 required arguments and 4 optional arguments:
- *--fi* - **required** argument for path to first cell line gene expressions table (**.csv** format only)
- *--si* - **required** argument for path to second cell line gene expressions table (**.csv** format only)
- *--out* - **required** argument for path for the output file (will be written in **.csv** format)
- *--need_ci* - **optional** argument, if you write it in your command, the result of confidence intervals intersection in terms True/False will be added to the output
- *--adj* - **optional** argument, if you write it in your command, then your p-values (z-test) will be **multiple corrected** with default method bonferroni and default alpha 0.05 using **multipletests** from **statsmodels.ststs.multitest** library
- *--adj_method* - **optional** argument to specify any other multiple correction method
- *--adj_alpha* - **optional** argument to specify any other alpha value for multiple correction method

## List of correction methods that could be chosen for --adj_method argument ##

- bonferroni
- sidak 
- holm-sidak
- holm
- simes-hochberg
- hommel
- fdr_bh 
- fdr_by
- fdr_tsbh
- fdr_tsbky

## Other functions of the script ##

The `diffexpr` function uses the following functions: 
- `diffexpr_ci` - takes two tables with gene expressions and gene names list and returns a list of boolean values that are answers to a question "whether the confidence intervals of gene expression between two cell lines do not intersect and we can assume the difference to be significant?"
- `diffexpr_ztest` - takes two tables with gene expressions and gene names list and returns a tuple of two lists: list of boolean values if the difference is significant and list of p-values
- `diff_means` - takes two tables with gene expressions and gene names list and returns a list of differences between mean expression of each gene between second and first tables
- `expression_means` - takes one table and gene names list and returns a list of means for expression of each gene in a table

## Output file ##

The output file is a table in **.csv** format that can contain:
- gene names in alhpabetic order (always)
- mean expression of each gene in first table (always)
- mean expression of each gene in second table (always)
- difference between mean gene expresison in the second cell line and first cell line (always)
- results of ci-test (only if chosen; True if gene expression differs significantly between cell lines, False otherwise)
- z-test p-value (always, not corrected for multiple comparisons)
- z-test adjusted p-value (only if adjustment chosen)
- test outcome (True if gene expression differs significantly between cell lines, False otherwise; calculated according non adjusted z-test p-value if no adjustment chosen, calculated according adjusted z-test p-value otherwise)

### To run the script correctly please use **python 3** ###
### This script uses several libraries, so **please install** their versions from the *requirements.txt* file to your virtual environment ###

## Script usage example ##

I recommend you to run the script from shell:

` python3 diffexpress.py --fi {path/to/the/first/table.csv} --si {path/to/the/second/table.csv} --out {name_of_output_table} --ci --adj --adj_method sidak --adj_alpha 0.001`

Remember that only *--fi*, *--si* and *--out* are required arguments, other are optional.
