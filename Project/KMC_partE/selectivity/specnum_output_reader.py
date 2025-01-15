import pandas as pd

def sum_nonzero_values(file_path, column_name):
    try:
        # Read the file with space as delimiter
        df = pd.read_csv(file_path, delimiter=r"\s+", engine="python")
    except Exception as e:
        print(f"Error reading the file: {e}")
        return

    # Strip spaces from column names
    df.columns = df.columns.str.strip()

    # Check if the column exists
    if column_name not in df.columns:
        print(f"Column '{column_name}' not found in the file.")
        print(f"Available columns: {df.columns.tolist()}")
        return

    try:
        # Convert the column to numeric, coercing errors
        df[column_name] = pd.to_numeric(df[column_name], errors='coerce')

        # Filter out non-zero values and calculate their sum
        sum_nonzero = df.loc[df[column_name] != 0, column_name].sum()
        print(f"Sum of non-zero values in column '{column_name}': {sum_nonzero}")
        return sum_nonzero
    except Exception as e:
        print(f"Error processing column '{column_name}': {e}")
        return

# Example usage
file_path = "/Users/lotburgstra/Desktop/TCC_SM_/Project/KMC/selectivity/specnum_output.txt"
column_name = "CH4*"
sum_nonzero_values(file_path, column_name)
