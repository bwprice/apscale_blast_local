#!/bin/bash
#SBATCH --partition=hour
#SBATCH --output=%j_apscale_out.out
#SBATCH --error=%j_apscale_err.err
#SBATCH --mem=40G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=name@email.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=apscale_blast

# Load required modules
module load blast/2.9.0  # Load the BLAST module

# Set variables (modify these for your specific run)
DATABASE="/gpfs/nhmfsa/bulk/share/data/mbl/share/workspaces/users/benjp/software/apscale_blast/MIDORI/"  # Path to your BLAST database
QUERY_FASTA="/gpfs/nhmfsa/bulk/share/data/mbl/share/workspaces/users/benjp/software/apscale_blast/big_test.fasta"  # Path to your query sequences
OUTPUT_DIR="/gpfs/nhmfsa/bulk/share/data/mbl/share/workspaces/users/benjp/software/apscale_blast/output/"  # Where to save results
N_CORES=16  # Should match --cpus-per-task above

# Optional parameters
TASK="blastn"  # Options: blastn, megablast, dc-megablast
SUBSET_SIZE=100  # Number of sequences per subset
MAX_TARGET_SEQS=20  # Number of hits to retain
THRESHOLDS="97,95,90,87,85"  # Taxonomy filter thresholds

echo "Starting apscale_blast job at $(date)"
echo "Database: $DATABASE"
echo "Query: $QUERY_FASTA"
echo "Output: $OUTPUT_DIR"
echo "Cores: $N_CORES"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run apscale_blast with your modified version
apscale_blast \
    -db "$DATABASE" \
    -q "$QUERY_FASTA" \
    -o "$OUTPUT_DIR" \
    -n_cores "$N_CORES"
    -task "$TASK" \
    -subset_size "$SUBSET_SIZE" \
    -max_target_seqs "$MAX_TARGET_SEQS" \
    -thresholds "$THRESHOLDS"

# Check if the job completed successfully
if [ $? -eq 0 ]; then
    echo "apscale_blast completed successfully at $(date)"
    echo "Results saved to: $OUTPUT_DIR"
else
    echo "apscale_blast failed with exit code $?" >&2
    exit 1
fi

echo "Job complete!"
