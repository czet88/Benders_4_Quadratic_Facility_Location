import pandas as pd
import sys
import time
import os
from glob import glob


def create_latex_table(df, col_list=None):
    # Function that reads the dataframe and a list of ordered columns
    # to be included and outputs the string to be pasted in a latex document
    if not df.empty:
        if col_list:
            cols = col_list
        else:
            cols = df.columns.tolist()
        out_string = df[cols].to_latex(index=False)
    else:
        out_string = "Empty dataframe"
    return out_string


def from_csvs_to_dataframe(file_prefix):
    # Function that reads all files in the current folder beginning with file_prefix
    # and concatenates them all into one dataframe.
    print("Searching for files beginning with " + file_prefix + ".")
    path_pf = f'./{file_prefix}*.csv'
    file_list = glob(path_pf)
    list_of_dfs = []
    for file in file_list:
        list_of_dfs.append(pd.read_csv(file, delimiter=";"))
    full_df = pd.concat(list_of_dfs)
    return full_df


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

    # List of conditional statements
    uncap_cond = df_full['Cap'] == 0
    cap_cond = df_full['Cap'] == 1
    fix_charge_cond = df_full['hybrid'] == 1
    p_med_cond = df_full['hybrid'] == 0

    # Uncapacitated dataframes
    df_uncap_fix_charge = df_full[uncap_cond & fix_charge_cond]
    df_uncap_p_med = df_full[uncap_cond & p_med_cond]

    # Capacitated dataframes
    df_cap_fix_charge = df_full[cap_cond & fix_charge_cond]
    df_cap_p_med = df_full[cap_cond & p_med_cond]

    # Create the Latex tables
    with open(outfile_name, "w") as text_file:
        print(f"{create_latex_table(df_uncap_fix_charge, ['Vers', 'instance'])}", file=text_file)
        print(f"{create_latex_table(df_uncap_p_med, ['Vers', 'instance'])}", file=text_file)
        print(f"{create_latex_table(df_cap_fix_charge, ['Vers', 'instance'])}", file=text_file)
        print(f"{create_latex_table(df_cap_p_med, ['Vers', 'instance'])}", file=text_file)


if __name__ == "__main__":
    main()
