# Codon Usage Analyzer and Optimizer

## What does this project do?

This project is a tool that analyzes the codon usage in a DNA sequence and suggests an optimized version based on the preferred codon usage of a target organism (e.g., *E. coli* or yeast). It helps improve gene expression by recommending codon changes that do not alter the protein but make the gene easier for the organism to translate efficiently.

The tool also provides useful information such as GC content and highlights rare codons that might reduce protein production efficiency. It visualizes codon usage with simple charts to help the user understand the sequence’s properties.

## Input and Output

- **Input:**  
  - A DNA sequence file (FASTA or plain text format) and a target organism's NCBI Taxonomy ID for codon optimization.

- **Output:**   
  - An optimized DNA sequence with codons adjusted to the target organism’s preferences (FASTA format).  
  - Summary statistics including GC content, rare codon warnings, and codon frequencies (text file).  
  - Visual charts (PNG) showing codon usage before and after optimization.  

## Technical details

- Written in Python, using libraries such as Biopython for sequence parsing and Matplotlib for visualization.
- Runs via command line interface (CLI).
- Input files are read locally; outputs saved to a project folder.
- Dependencies can be installed via `pip install -r requirements.txt`.

## How to download, install, and run

1. Clone the repository:
   ```bash
   git clone https://github.com/leehalevi/Codon-Usage-Analyzer-Optimizer.git
   cd Codon-Usage-Analyzer-Optimizer
   ```
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Run the program with command-line arguments:
   ```bash
   python codon_optimizer.py --input mygene.fasta --taxid 316407
   ```
   All output files will be saved in the `results` directory.

   - `--input`: Input DNA sequence file (FASTA or plain text)
   - `--taxid`: Target organism NCBI Taxonomy ID (see list below)
   - `--output`: Output file for optimized sequence (default: optimized_sequence.fasta). This will be saved in the `results` directory.
   - `--stats`: Output file for statistics (default: sequence_stats.txt). This will be saved in the `results` directory.

The following NCBI Taxonomy IDs are supported out-of-the-box:

- Bacillus subtilis: 1423
- Caenorhabditis elegans: 6239
- Drosophila melanogaster: 7227
- Escherichia coli: 316407
- Gallus gallus (chicken): 9031

