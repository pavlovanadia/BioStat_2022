import pandas as pd
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest
import argparse


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


def diffexpr(
    first_cell_type_expressions_path,
    second_cell_type_expressions_path,
    save_results_table
    ):
    
    table_1 = pd.read_csv(first_cell_type_expressions_path, index_col=0) # upload first table
    table_2 = pd.read_csv(second_cell_type_expressions_path, index_col=0) # upload second table

    # gene name has to be in both tables and the cell type has not to be considered as a gene name
    gene_names = list(set(table_1.columns).intersection(set(table_2.columns)).difference(set(["Cell_type"])))

    # list with gene names and ci results
    ci_results = diffexpr_ci(table_1, table_2, gene_names)

    # list with ztest results and ztest p values
    z_results, z_p_vals = diffexpr_ztest(table_1, table_2, gene_names)

    # list with difference between means of expressions
    mean_diff = diff_means(table_1, table_2, gene_names)

    # creating of a dataframe
    results = {
        "gene_name": gene_names,
        "ci_test_results": ci_results,
        "z_test_results": z_results,
        "z_test_p_values": z_p_vals,
        "mean_diff": mean_diff
    } 

    results = pd.DataFrame(results)
    
    results.to_csv(f"{save_results_table}.csv", index=False)


parser = argparse.ArgumentParser(description='Differential expression analysis') # created parser
parser.add_argument('first_inp', type=str, help='Path to first (control) .csv file') # argument for first table
parser.add_argument('second_inp', type=str, help='Path to second (experiment) .csv file') # argument for second table
parser.add_argument('output_name', type=str, help='Results table name') # argument for output
args = parser.parse_args()

# finally let the script do the work
diffexpr(
    first_cell_type_expressions_path=args.first_inp,
    second_cell_type_expressions_path=args.second_inp,
    save_results_table=args.output_name
)
