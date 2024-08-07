import sys
import pandas as pd

def extract_columns(matrix_file, column_names, output_file):
    matrix_df = pd.read_csv(matrix_file, sep = "\t", dtype = "str", index_col=0)
    extracted_columns = matrix_df[column_names]
    extracted_columns.to_csv(output_file, sep='\t')

def exclude_columns(matrix_file, exclude_columns, output_file):
    matrix_df = pd.read_csv(matrix_file, sep = "\t", dtype = "str", index_col=0)
    columns_to_extract = [col for col in matrix_df.columns if col not in exclude_columns]
    extracted_columns = matrix_df[columns_to_extract]
    extracted_columns.to_csv(output_file, sep='\t')

if __name__ == "__main__":
    matrix_file = sys.argv[1] 
    column_names = sys.argv[2].split(",")
    output_file = sys.argv[3]
    if len(sys.argv) > 4 and sys.argv[4] == "-v":
        exclude_columns(matrix_file, column_names, output_file)
        print("Columns",  sys.argv[2], " EXCLUDED and the rest is saved to ", output_file)
    else:
        extract_columns(matrix_file, column_names, output_file)
        print("Columns ",  sys.argv[2], " extracted and saved to ", output_file)


