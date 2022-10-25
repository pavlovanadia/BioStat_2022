# Differential gene expression #

This script `diffexpress.py` has several functions that calculate the difference between gene expressions of two cell lines.

**The main function of the script** is `diffexpr` that takes 3 arguments:
- first_cell_type_expressions_path - path to first cell line gene expressions table (**.csv** format only)
- second_cell_type_expressions_path - path to second cell line gene expressions table (**.csv** format only)
- save_results_table - path for the output file (**.csv** format)

The `diffexpr` function uses the following functions: 
- `diffexpr_ci`
- `diffexpr_ztest`
- `diff_means`

The output file is a table in **.csv** format that contains:
- gene name
- results of ci-test (True if gene expression differs significantly between cell lines, False otherwise)
- results of z-test (True if ztest expression differs significantly between cell lines, False otherwise)
- z-test p-value
- difference between mean gene expresison in the second cell line and first cell line

### To run the script correctly please use **python 3** ###
### This script uses numpy and pandas libraries, so **please install** their versions from the *requirements.txt* file to your virtual environment ###