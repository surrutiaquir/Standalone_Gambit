import pandas as pd
import numpy as np
import sys

def process_csv(filename):
    # Load the CSV file and select the first two columns
    df = pd.read_csv(filename, header=None, usecols=[0, 1], names=['logM', 'logU'])
    
    # Calculate 10^logM and 10^logU
    df['M'] = 10 ** df['logM']
    df['U'] = 10 ** df['logU']

    # Check if the 'M' column is strictly increasing
    if (df['M'].diff().dropna() > 0).all():
        print("The 'M' column is strictly increasing.")
    else:
        print("The 'M' column is not strictly increasing. Making it strictly increasing by sorting and removing duplicates.")
        
        # Sort the DataFrame by 'M' in ascending order
        df = df.sort_values(by='M').reset_index(drop=True)
        
        # Drop duplicate 'M' values, keeping the first occurrence
        df = df.drop_duplicates(subset='M', keep='first').reset_index(drop=True)

    # Create a new filename for the processed CSV
    new_filename = filename.rsplit('.', 1)[0] + '_proc.csv'

    # Save the new CSV file
    df[['M', 'U']].to_csv(new_filename, index=False, header=False)

    # Print the min and max values of 10^X
    min_value = df['M'].min()
    max_value = df['M'].max()
    
    print(f"Min value of M_N (GeV): {min_value}")
    print(f"Max value of M_N (GeV): {max_value}")

# Example usage
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
    else:
        process_csv(sys.argv[1])
