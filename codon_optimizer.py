from Bio import SeqIO
from Bio.Data import CodonTable
import python_codon_tables as pct
import argparse
import matplotlib.pyplot as plt
import os

#Available tables: ['b_subtilis_1423', 'c_elegans_6239', 'd_melanogaster_7227', 'e_coli_316407', 'g_gallus_9031', 'h_sapiens_9606', 'm_musculus_10090', 'm_musculus_domesticus_10092', 's_cerevisiae_4932']
# Also from a taxID

def get_codon_usage_table(taxid):
    """
    Returns a codon usage frequency dictionary for the given taxid using python_codon_tables.
    Output: dict {codon (str): frequency (float)}
    """
    table = pct.get_codons_table(taxid)
    codon_usage = {}
    # The table is a nested dict: {amino_acid: {codon: freq, ...}}
    for aa_codons in table.values():
        for codon, freq in aa_codons.items():
            codon_usage[codon.upper()] = freq
    return codon_usage

def read_sequence(input_file):
    # Try to parse as FASTA
    try:
        record = next(SeqIO.parse(input_file,"fasta"))
        return str(record.seq)
    except Exception:
        # Fallback: read as plain text
        with open(input_file) as f:
            return f.read().replace("\n","").strip() #should this be error regulated as well?

def count_codon_frequencies(sequence):
    codon_counts = {}
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3].upper()
        if len(codon) == 3:
            codon_counts[codon] = codon_counts.get(codon, 0) + 1
    return codon_counts

def calculate_gc_content(sequence):
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100 if sequence else 0

def identify_rare_codons(sequence, codon_usage_table, threshold=0.1):
    """
    Identify rare codons in the sequence based on the codon usage table.
    Args:
        sequence (str): DNA sequence.
        codon_usage_table (dict): Codon usage frequencies (codon -> frequency).
        threshold (float): Codons with usage below this fraction are considered rare.
    Returns:
        list: List of rare codons found in the sequence.
    """
    codon_freqs = count_codon_frequencies(sequence)
    rare_codons = []
    for codon in codon_freqs:
        usage = codon_usage_table.get(codon, 0)
        if usage < threshold:
            rare_codons.append(codon)
    return rare_codons

def optimize_sequence(sequence, codon_usage_table):
    """
    Returns a codon-optimized DNA sequence for the target organism.
    For each amino acid, replaces the codon with the most frequent codon in the codon usage table.
    """
    from Bio.Seq import Seq
    from Bio.Data import CodonTable
    # Build a mapping: amino acid -> preferred codon
    aa_to_best_codon = {}
    for codon, freq in codon_usage_table.items():
        aa = str(Seq(codon).translate(table=11))  # Standard table
        if aa not in aa_to_best_codon or freq > codon_usage_table[aa_to_best_codon[aa]]:
            aa_to_best_codon[aa] = codon
    # Translate input sequence to amino acids
    seq_obj = Seq(sequence)
    aa_seq = str(seq_obj.translate(to_stop=False, table=11))
    # Build optimized sequence
    optimized_codons = []
    for aa in aa_seq:
        if aa == '*':  # Stop codon
            # Find the most frequent stop codon
            stop_codons = [c for c in codon_usage_table if str(Seq(c).translate(table=11)) == '*']
            best_stop = max(stop_codons, key=lambda c: codon_usage_table[c])
            optimized_codons.append(best_stop)
        else:
            optimized_codons.append(aa_to_best_codon.get(aa, 'NNN'))
    return ''.join(optimized_codons)

def plot_codon_usage(codon_freqs, title, filename):
    # Ensure results directory exists
    results_dir = os.path.dirname(filename)
    if results_dir and not os.path.exists(results_dir):
        os.makedirs(results_dir)
    codons = list(codon_freqs.keys())
    counts = [codon_freqs[c] for c in codons]
    plt.figure(figsize=(18, 5))
    plt.bar(codons, counts)
    plt.title(title)
    plt.xlabel('Codon')
    plt.ylabel('Count')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Codon Usage Analyzer and Optimizer")
    parser.add_argument('--input', required=True, help='Input DNA sequence file (FASTA or plain text)')
    parser.add_argument('--taxid', required=True, type=int, help='Target organism NCBI Taxonomy ID')
    parser.add_argument('--output', default='optimized_sequence.fasta', help='Output file for optimized sequence (FASTA)')
    parser.add_argument('--stats', default='sequence_stats.txt', help='Output file for statistics')
    args = parser.parse_args()

    input_file = args.input
    taxid = args.taxid

    results_dir = "results"
    os.makedirs(results_dir, exist_ok=True)

    output_file = os.path.join(results_dir, os.path.basename(args.output))
    stats_file = os.path.join(results_dir, os.path.basename(args.stats))

    sequence = read_sequence(input_file)
    codon_freqs = count_codon_frequencies(sequence)
    gc_content = calculate_gc_content(sequence)

    print("Codon frequencies:", codon_freqs)
    print(f"GC content: {gc_content:.2f}%")

    codon_usage_table = get_codon_usage_table(taxid)
    if not codon_usage_table:
        print(f"Could not retrieve codon usage table for taxid: {taxid}")
        return
    print(f"Codon usage table: {codon_usage_table}")

    rare_codons = identify_rare_codons(sequence, codon_usage_table)
    print("Rare codons in the organism that appear in this sequence:", rare_codons)

    optimized_seq = optimize_sequence(sequence, codon_usage_table)
    print(f"Optimized sequence (saved into {output_file}):", optimized_seq)
    with open(output_file, "w") as f:
        f.write(">optimized_sequence\n")
        for i in range(0, len(optimized_seq), 60):
            f.write(optimized_seq[i:i+60] + "\n")

    # Save statistics
    with open(stats_file, "w") as f:
        f.write(f"Input file: {input_file}\n")
        f.write(f"TaxID: {taxid}\n")
        f.write(f"GC content: {gc_content:.2f}%\n")
        f.write(f"Rare codons: {rare_codons}\n")
        f.write(f"Codon frequencies: {codon_freqs}\n")

    # Save gene_stats and optimized_gene in results directory
    with open(os.path.join(results_dir, "gene_stats.txt"), "w") as f:
        f.write(f"GC content: {gc_content:.2f}%\n")
        f.write(f"Rare codons: {rare_codons}\n")
        f.write(f"Codon frequencies: {codon_freqs}\n")
    with open(os.path.join(results_dir, "optimized_gene.fasta"), "w") as f:
        f.write(">optimized_gene\n")
        for i in range(0, len(optimized_seq), 60):
            f.write(optimized_seq[i:i+60] + "\n")

    # Visualization
    plot_codon_usage(codon_freqs, "Original Codon Usage", os.path.join(results_dir, "original_codon_usage.png"))
    optimized_freqs = count_codon_frequencies(optimized_seq)
    plot_codon_usage(optimized_freqs, "Optimized Codon Usage", os.path.join(results_dir, "optimized_codon_usage.png"))
    print(f"Codon usage plots saved as {os.path.join(results_dir, 'original_codon_usage.png')} and {os.path.join(results_dir, 'optimized_codon_usage.png')}")

if __name__ == '__main__':
    main()
