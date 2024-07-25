#!/bin/bash

import csv
import os

# define the path to the input CSV file
input_csv_path = 'Path/to/updated_mirnas_of.interest.csv'

# define and create the directory to store the output files
output_dir = 'Desired/path/to/mirnas_of_interest_482''
os.makedirs(output_dir, exist_ok=True)

# read the input CSV file 
# content flagged with '[!]' was co-piloted with ChatGPT, an AI language model by OpenAI.
with open(input_csv_path, newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        content = row[0] # [!]
        filename = f"{content}.csv" # [!]
        file_path = os.path.join(output_dir, filename) #[!]

        # write the content to the new file [!]
        with open(file_path, 'w', newline='') as newfile:
            writer = csv.writer(newfile)
            writer.writerow([content])

print(f'Written .csv file can be found in {output_dir}')
