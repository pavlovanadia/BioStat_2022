import pandas as pd
import numpy as np
import argparse
import sys

import scipy.stats as st
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.weightstats import ztest


def check_intervals_intersect(first_ci, second_ci):
    """If True, then the intervals do intersect and therefore the difference is not significant
    If False, then the intervals do not intersect and therefore the difference is significant"""
    return max(first_ci[0], second_ci[0]) < min(first_ci[1], second_ci[1])


def diffexpr_ci(first_table, second_table, gene_names):
    """function for a tool that calculates the difference between gene expressions from two tables using confidence intervals"""

    ci_test_results = list() # results list

    for gene in gene_names:
        ci_first_table = st.t.interval(
            alpha=0.95, # confidence interval for gene expression in cells from fisrt table
            df=len(first_table[gene]) - 1,
            loc=np.mean(first_table[gene]),
            scale=st.sem(first_table[gene]))

        ci_second_table = st.t.interval(
            alpha=0.95, # confidence interval for gene expression in cells from second table
            df=len(second_table[gene]) - 1, 
            loc=np.mean(second_table[gene]), 
            scale=st.sem(second_table[gene]))

        # if intervals do not intersect we consider the difference significant
        ci_test_results.append(not (check_intervals_intersect(ci_first_table, ci_second_table))) 
    return ci_test_results # two lists so that gene name would correspond to ci test result


def diffexpr_ztest(table_1, table_2, gene_names):
    """function for a tool to calculate differences between gene expressions usinf ztest"""

    z_test_values = list() # list for p_values
    z_test_results = list() # list for results in terms True/False 

    for gene in gene_names:
        z_result = ztest(
            table_1[gene], 
            table_2[gene]
            )
        p_val = z_result[1]
        z_test_values.append(p_val)
        z_test_results.append(p_val < 0.05)
    
    return z_test_results, z_test_values


def diff_means(table_1, table_2, gene_names):
    """calculate the difference benween means of each gene expressions between two tables"""

    difference = list()

    for gene in gene_names:
        difference.append(np.mean(table_2[gene]) - np.mean(table_1[gene]))
    return difference


def expression_means(table, gene_names):
    """calculate the means of each gene expression for one table"""

    expr_means = list()

    for gene in gene_names:
        expr_means.append(np.mean(table[gene]))

    return expr_means


def diffexpr(
    first_cell_type_expressions_path,
    second_cell_type_expressions_path,
    save_results_table,
    need_ci,
    need_adj,
    adj_method,
    adj_alpha
    ):
    
    table_1 = pd.read_csv(first_cell_type_expressions_path, index_col=0) # upload first table
    table_2 = pd.read_csv(second_cell_type_expressions_path, index_col=0) # upload second table

    # gene name has to be in both tables and the cell type has not to be considered as a gene name
    gene_names = sorted(list(set(table_1.columns).intersection(set(table_2.columns)).difference(set(["Cell_type"]))))

    # genes expression means for cell type from first table
    table_1_expressions = expression_means(table_1, gene_names)

    # genes expression means for cell type from second table
    table_2_expressions = expression_means(table_2, gene_names)

    # list with difference between means of expressions
    mean_diff = diff_means(table_1, table_2, gene_names)

    ci_results = None
    if need_ci: # if user wants to see the ci intersection (True/False) in output
        
        # list with gene names and ci results
        ci_results = diffexpr_ci(table_1, table_2, gene_names)

    # list with ztest p values
    outcome, z_p_vals = diffexpr_ztest(table_1, table_2, gene_names)

    # if p-values multiple comparisons correction is needed 
    z_p_vals_adjusted = None
    if need_adj:
        outcome, z_p_vals_adjusted, _, _ = multipletests(z_p_vals, adj_alpha, adj_method)

    # creating of a dataframe
    results = {
        "gene_name": gene_names,
        "first_cell_type_expression_means": table_1_expressions,
        "second_cell_type_expressions_means": table_2_expressions,
        "mean_difference": mean_diff,
        "ci_test_results": ci_results,
        "z_test_p_values": z_p_vals,
        f"p_values_adjusted_method_{adj_method}_alpha_{adj_alpha}": z_p_vals_adjusted,
        "test_outcome": outcome
    }

    results = pd.DataFrame(results)
    results.dropna(how='all', axis=1, inplace=True) # if ci and multiple test correcrion not chosen, these columns will be deleted
    
    results.to_csv(f"{save_results_table}.csv", index=False)


""" List of multiple comparisons p-values adjustment methods """
FDR_METHODS_ALLOWED = [
    "bonferroni", 
    "sidak", 
    "holm-sidak", 
    "holm", 
    "simes-hochberg", 
    "hommel",
    "fdr_bh", 
    "fdr_by",
    "fdr_tsbh",
    "fdr_tsbky"
]

# arguments parsing from command line
parser = argparse.ArgumentParser(description='Differential expression analysis') # created parser
parser.add_argument('--fi', type=str, help='Path to first (control) .csv file') # argument for first table
parser.add_argument('--si', type=str, help='Path to second (experiment) .csv file') # argument for second table
parser.add_argument('--out', type=str, help='Results table name') # argument for output
parser.add_argument('--ci', action='store_true', dest='ci', help='Do you want to calculate if ci intersect') # flag for ci intersect
parser.add_argument('--adj', action='store_true', dest='adj', help='Do you want to correct p-velues for multiple comparisons?')
parser.add_argument('--adj_method', default="bonferroni", choices=FDR_METHODS_ALLOWED, help='Type p-value correction method for multiple comparisons')
parser.add_argument('--adj_alpha', default=0.05, type=float, help='alpha for multiple comparisons')

if len(sys.argv) < 2: # autimatically will print help if no arguments provided
    parser.print_help()
    sys.exit(0)
args = parser.parse_args()


# finally let the script do the work
diffexpr(
    first_cell_type_expressions_path=args.fi,
    second_cell_type_expressions_path=args.si,
    save_results_table=args.out,
    need_ci=args.ci,
    need_adj=args.adj,
    adj_method=args.adj_method,
    adj_alpha=args.adj_alpha
)
