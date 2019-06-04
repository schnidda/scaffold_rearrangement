# scaffold_rearrangement
scaffold rearrangement tool
v1.0
2019 06 03
written by Fabian Schneider

This script uses the Autodesk nano design package.
A modified version comes with this script, as I implemented some of the DIetzlab scaffold sequences manually to its library.
They can be found in nanodesign\nanodesign\converters\dna_sequence_data.py
If new sequences are needed or wished, paste them there.

Afterwards they can be used as template for new sequences.


To run the script:
python 2.7 is needed
have a .json file and a .txt file with provided staple sequences in the same folder as the script
start the script in your terminal using: python scaffold_rearrange_v1.py
The script asks you for .json input: "your_file.json“
The script asks you for the scaffold sequence you want to use as template: "pCS2"
The script asks you for provided staples: "your_staples.txt"

After successful run, you will find a file called created_scaffold.txt
