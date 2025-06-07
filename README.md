# Codon Usage Analyzer and Optimizer

## What does this project do?

This project is a tool that analyzes the codon usage in a DNA sequence and suggests an optimized version based on the preferred codon usage of a target organism (e.g., *E. coli* or yeast). It helps improve gene expression by recommending codon changes that do not alter the protein but make the gene easier for the organism to translate efficiently.

The tool also provides useful information such as GC content and highlights rare codons that might reduce protein production efficiency. It visualizes codon usage with simple charts to help the user understand the sequence’s properties.

## Input and Output

- **Input:**  
  - A DNA sequence file (FASTA or plain text format) and a choice of target organism for codon optimization.

- **Output:**  
  - An optimized DNA sequence with codons adjusted to the target organism’s preferences.  
  - Summary statistics including GC content and rare codon warnings.  
  - Visual charts showing codon usage before and after optimization.  

## Technical details

- Written in Python, using libraries such as Biopython for sequence parsing and Matplotlib or Plotly for visualization.  
- Runs via command line interface (CLI) or simple graphical interface (optional).  
- Input files are read locally; outputs saved to a project folder.  
- Dependencies can be installed via `pip install -r requirements.txt`.

## How to download, install, and run

1. Clone the repository:  
   ```bash
   git clone https://github.com/leehalevi/Codon-Usage-Analyzer-Optimizer.git
   cd Codon-Usage-Analyzer-Optimizer

2. Install dependencies:
   ```bash
   pip install -r requirements.txt

3. Run the program with an example command:
   ```bash
   python codon_optimizer.py --input mygene.fasta --organism ecoli --output optimized_gene.fasta

---

*This project was created as part of the course Basic Programming Skills in Python at the Weizmann Institute of Science.*


