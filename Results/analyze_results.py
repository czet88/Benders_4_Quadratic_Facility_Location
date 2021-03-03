import pandas as pd
import numpy as np
import re
import sys
import time
import os
from glob import glob


def create_latex_table(df, col_list=None, caption_string=None):
    # Function that reads the dataframe and a list of ordered columns
    # to be included and outputs the string to be pasted in a latex document
    if not df.empty:
        if col_list:
            cols = col_list
        else:
            cols = df.columns.tolist()
        # Add caption if there is one
        if caption_string:
            out_string = df[cols].to_latex(index=True, index_names=True, multicolumn_format='c', caption=caption_string)
        else:
            out_string = df[cols].to_latex(index=True, index_names=True, multicolumn_format='c')
    else:
        out_string = "Empty dataframe"
    return out_string


def name_heur(row):
    # Create a column with the name of the heuristic used
    if row['heurParam'] == 1:
        return 'MH'
    if row['heurParam'] == 2:
        return 'MH_LS'
    if row['heurParam'] == 3:
        return 'MH_LS_ILS'
    return 'Other'


def from_csvs_to_dataframe(file_prefix):
    # Function that reads all files in the current folder beginning with file_prefix
    # and concatenates them all into one dataframe.
    print("Searching for files beginning with " + file_prefix + ".")
    path_pf = f'./{file_prefix}*'
    file_list = glob(path_pf)
    list_of_dfs = []
    for file in file_list:
        list_of_dfs.append(pd.read_csv(file, delimiter=";"))
    full_df = pd.concat(list_of_dfs)

    # Rename some of the columns
    full_df['Heur_used'] = full_df.apply(lambda row: name_heur(row), axis=1)
    full_df.rename(columns={'instance': 'Instance'}, inplace=True)

    # Add columns for additional reporting
    full_df['Rel_Heur_Dev'] = 100*(full_df['UBPre']-full_df['UB']) / full_df['UB']
    full_df['Rel_Time_Root'] = 100 * (full_df['CPU_Pre']) / full_df['CPU_all']
    return full_df


def split_dataframe_by_case(df_full):
    # List of conditional statements
    uncap_cond = df_full['Cap'] == 0
    cap_cond = df_full['Cap'] == 1
    fix_charge_cond = df_full['hybrid'] == 1
    p_med_cond = df_full['hybrid'] == 0

    df_dict = {}

    # Uncapacitated dataframes
    df_uncap_fix_charge = df_full[uncap_cond & fix_charge_cond]
    df_uncap_p_med = df_full[uncap_cond & p_med_cond]

    # Capacitated dataframes
    df_cap_fix_charge = df_full[cap_cond & fix_charge_cond]
    df_cap_p_med = df_full[cap_cond & p_med_cond]

    # Add to the dictionary only if it is not empty
    if not df_uncap_fix_charge.empty:
        df_dict["uncap_fixed_cost"] = [df_uncap_fix_charge, "Results for the uncapacitated fixed charge problems"]
    if not df_uncap_p_med.empty:
        df_dict["uncap_p_med"] = [df_uncap_p_med, "Results for the uncapacitated p-median problems"]
    if not df_cap_fix_charge.empty:
        df_dict["cap_fixed_cost"] = [df_cap_fix_charge, "Results for the capacitated fixed charge problems"]
    if not df_cap_p_med.empty:
        df_dict["cap_p_med"] = [df_cap_p_med, "Results for the capacitated p-median problems"]

    return df_dict


def add_report_table(df_dict, report_vals, index_vals, col_vals):
    for ele in df_dict:
        # If we are using the p-median problem, then we need to include the value of p
        if re.search('p_med', ele):
            index_app = ['p']
        else:
            index_app = []
        index_use = index_vals + index_app
        df_rep = pd.pivot_table(df_dict[ele][0], values=report_vals, index=index_use, columns=col_vals, aggfunc=np.mean,
                                margins=True, margins_name='Testbed average').round(2)
        new_order = ['CPU_all','Rel_Heur_Dev','Num_fixed','Rel_Time_Root','BBnodes']
        df_rep = df_rep.reindex(new_order, axis=1)
        # df_rep = df_rep.iloc[:, :-1]
        # df_rep = df_rep.drop('Testbed average', axis=1, level=1)
        df_dict[ele].append(df_rep)

    return df_dict


def create_latex_tables(df_dict, df_position, caption_position, outfile_name):
    # Function that receives a dictionary of dataframes and captions and prints the corresponding latex tables to
    # outfile_name
    with open(outfile_name, "w") as text_file:
        for ele in df_dict:
            print(f"{create_latex_table(df_dict[ele][df_position], caption_string=df_dict[ele][caption_position])}",
                  file=text_file)
    return None


def main():
    # Setup
    output_folder_name = "Latex"
    file_prefix = "R_"
    time_str = time.strftime("%Y%m%d_%H%M")
    outfile_name = f"./{output_folder_name}/Ltx_{time_str}.txt"

    # Ensure output folder exists
    if not os.path.exists(output_folder_name):
        os.makedirs(output_folder_name)

    # Read arguments and adjust if necessary
    if len(sys.argv) >= 2:
        file_prefix = sys.argv[1]

    # Create one concatenated dataframe for all tables read
    df_full = from_csvs_to_dataframe(file_prefix)

    # Split the dataframe into each of the cases
    df_dict = split_dataframe_by_case(df_full)

    # Define the columns, indices and values to be reported from the dataframes in the dictionaries
    report_vals = ['CPU_all','Rel_Heur_Dev','Num_fixed','Rel_Time_Root','BBnodes']
    index_vals = ['Instance']
    col_vals = None  # ['Heur_used']

    # Add the dataframes with the reports to the dataframe dictionary
    df_dict = add_report_table(df_dict, report_vals, index_vals, col_vals)

    # Create the Latex tables
    create_latex_tables(df_dict, 2, 1, outfile_name)

    print("Finished writing latex tables into " + str(outfile_name))


if __name__ == "__main__":
    main()
