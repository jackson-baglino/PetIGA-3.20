import os
import re

def rename_files_sequentially(directory):
    # List all files in the directory
    files = os.listdir(directory)

    # Regular expression to match the files with the original pattern
    pattern = re.compile(r'solV_(\d+)_contour\.vtp$')

    # Create a list of tuples (original_name, original_number)
    sorted_files = []
    for filename in files:
        match = pattern.match(filename)
        if match:
            original_number = int(match.group(1))
            sorted_files.append((filename, original_number))

    # Sort the list of tuples based on the original number
    sorted_files.sort(key=lambda x: x[1])

    # Loop through the sorted list and rename each file sequentially
    for i, (original_name, original_number) in enumerate(sorted_files):
        new_name = f"sol_contour_{i:03d}.vtp"
        os.rename(os.path.join(directory, original_name), os.path.join(directory, new_name))
        print(f"Renamed {original_name} to {new_name}")

# Example usage
rename_files_sequentially('./stlOut-copy/')