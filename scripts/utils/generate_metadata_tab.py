import os
import argparse
import pandas as pd
import ast
import sys

def get_unaligned_subfolders(config_file_path):
    try:
        # Read the content of the config file
        config_dict = {}
        exec(open(config_file_path).read(), config_dict)
        print(config_dict["unaligned"])
        if 'unaligned' in config_dict and isinstance(config_dict["unaligned"], list):
            unaligned_list = config_dict["unaligned"]

            # Extract subfolder names from the path in the list
            subfolder_names = []
            for path in config_dict["unaligned"]:
                subfolder_names = subfolder_names + [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))] 
            return subfolder_names
        else:
            print("Error: 'unaligned' key not found or not a list in the config file.")
            return None

    except Exception as e:
        print(f"Error: {e}")
        return None

def generate_rows(directory_path, template_file_path,sample_name, output_file_path):
    # Get subfolder names from the "libraries.csv" file
    subfolder_names = get_unaligned_subfolders(os.path.join(directory_path, "config.py"))

    if subfolder_names is None:
        return

    # Read header and template lines from the template file
    with open(template_file_path, "r") as template_file:
        lines = template_file.readlines()

        if len(lines) >= 2:
            header = lines[0]
            template_line = lines[1]
        else:
            print("Error: The template file must contain at least two lines (header and template).")
            return

    if sample_name not in template_line:
        print(f"\nError: The template file must contain {sample_name}\n")
        sys.exit(1)
    # Create the output file and write the header
    with open(output_file_path, "w") as output_file:
        output_file.write(header)

        # Write a row for each subfolder
        for subfolder in subfolder_names:
            row = template_line.replace(sample_name, subfolder)
            output_file.write(row)

    print(f"Output file '{output_file_path}' generated successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate rows based on a template for each subfolder found in the input folder.")
    parser.add_argument("input_folder", help="Path to the input folder containing files with the pattern '*._libraries.csv'.")
    parser.add_argument("template_file", help="Path to the file containing the header and template lines.")
    parser.add_argument("sample_name", help="Sample ")

    args = parser.parse_args()

    output_file_path = "output_file.txt"  # You can change the output file name if needed
    generate_rows(args.input_folder, args.template_file, args.sample_name, output_file_path)

