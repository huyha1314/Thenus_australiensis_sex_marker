import re
import random

# Seed for reproducibility
random.seed(42)

with open("merged.aln.treefile", "r") as f:
    newick = f.read()

# Find all occurrences of "):" and replace them with ")boot:"
# where boot is a random integer between 70 and 100
def replace_node(match):
    boot = random.randint(70, 100)
    # Let's write some high support (>=95), some medium support (80-94), and some low support (<80)
    return f"){boot}:"

newick_mocked = re.sub(r"\):", replace_node, newick)

with open("merged_mock_boot.treefile", "w") as f:
    f.write(newick_mocked)

print("Mock treefile created with bootstrap values.")
