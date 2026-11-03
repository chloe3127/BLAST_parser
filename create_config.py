import shutil
import sys
import os

def main():
    if len(sys.argv) < 2:
        print("error: please provide a gene name as an argument for new config file")
        sys.exit(1)

    gene_name = sys.argv[1]
    template_file = "configs/config_template.py"
    new_file = f"configs/config_{gene_name}.py"

    if os.path.exists(new_file):
        print(f"Warning: {new_file} already exists. Overwrite? (y/n)")
        choice = input().lower()
        if choice != "y":
            print("Cancelled")
            sys.exit(0)
    
    shutil.copy(template_file, new_file)
    print(f"Created new config file: {new_file}. Please update with appropriate values for {gene_name}")

if __name__ == "__main__":
    main()