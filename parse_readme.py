import re
import os

with open("README.md", "r") as f:
    content = f.read()

# Define the script sections
sections = [
    ("01_fastp", "1. Quailty control of raw reads with fastp"),
    ("02_genomescope", "2. Genome survey using kmer distribution"),
    ("03_genome_assembly", "3. Genome assembly"),
    ("04_gene_prediction", "4.Gene structure prediction"),
    ("05_functional_annotation", "5. Functional annotation"),
    ("06_ncrna_prediction", "6. Noncoding RNAs prediction"),
    ("07_sex_marker", "7. Finding candidate sex-specific"),
    ("08_genome_assessment", "9. Genome assessment"),
    ("09_phylogenetic_analysis", "10. Phylogenetic analysis")
]

os.makedirs("scripts", exist_ok=True)
os.makedirs("logs", exist_ok=True)

# Just splitting by ## or ### and extracting bash blocks
# Let's write a simple config
with open("config.sh", "w") as f:
    f.write('''#!/bin/bash
# Central configuration file for Thenus australiensis sex marker workflow

# Directories
export BASE_DIR="$PWD"
export DATA_DIR="$BASE_DIR/data"
export RESULT_DIR="$BASE_DIR/result"
export LOG_DIR="$BASE_DIR/logs"

# Ensure directories exist
mkdir -p "$DATA_DIR" "$RESULT_DIR" "$LOG_DIR"

# Compute resources
export THREADS=16
export MAX_MEMORY="80G"
''')

# We'll extract blocks based on headers
blocks = re.split(r'\n## ', '\n' + content)
script_idx = 1
for block in blocks:
    if not block.strip():
        continue
    title = block.split('\n')[0].strip()
    bash_matches = re.findall(r'```bash\n(.*?)\n```', block, re.DOTALL)
    if bash_matches:
        script_name = f"scripts/{script_idx:02d}_step.sh"
        # Try to make name more descriptive
        clean_title = re.sub(r'[^a-zA-Z0-9]+', '_', title.split('(')[0]).strip('_').lower()
        if len(clean_title) > 30:
            clean_title = clean_title[:30]
        script_name = f"scripts/{script_idx:02d}_{clean_title}.sh"
        
        with open(script_name, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("source config.sh\n\n")
            f.write('echo "Starting {}"\n'.format(title))
            for match in bash_matches:
                f.write(match + "\n")
        
        script_idx += 1

