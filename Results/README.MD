## Running the python script to produce csvs and Latex tables

Step 1: One must first execute the C code and copy the output files of the executions into the folder "Raw_Results"

Step 2: Open the terminal and execute the python script analyze_results.py as follows
        python analyze_results.py $FIRST_TWO_LETTERS_OF_OUTPUT_FILES
        
Step 3: The script will then produce a csv file for each problem variant(e.g. uncapacitated p median, etc.) and one
        file ending in .tex that contains in Latex format the tables of all variants found in the files that were
        analyzed.