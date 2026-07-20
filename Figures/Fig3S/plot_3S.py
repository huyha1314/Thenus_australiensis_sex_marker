import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sns
import numpy as np

# Set light, clean aesthetics
plt.style.use('default')
sns.set_context("notebook", rc={"font.size":12, "axes.titlesize":16, "axes.labelsize":14})
sns.set_style("white")

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Inter', 'Roboto', 'Arial']

# 1. LOAD DATA
df = pd.read_csv('/mnt/12T/lobster_project/script_new/test_sex_marker/2453/lobster_y_chromosome_candidates_separated.csv')
df = df[df['chrom'] == 'scaffold_2453'].copy()
df = df.sort_values('start')

# Smooth the depth data (rolling mean)
smooth_window = 3
df['wgs_m_smooth'] = df['wgs_m_median'].rolling(window=smooth_window, min_periods=1, center=True).mean()
df['wgs_f_smooth'] = df['wgs_f_median'].rolling(window=smooth_window, min_periods=1, center=True).mean()

# Smooth DArT-seq data
df['dart_m_smooth'] = df['dart_m_median'].rolling(window=smooth_window, min_periods=1, center=True).mean().fillna(0)
df['dart_f_smooth'] = df['dart_f_median'].rolling(window=smooth_window, min_periods=1, center=True).mean().fillna(0)

# 2. SETUP PLOT (2 PANELS)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
fig.patch.set_facecolor('white')
ax1.set_facecolor('white')
ax2.set_facecolor('white')

color_male_wgs = '#0072B2'     
color_female_wgs = '#D55E00'   
color_male_dart = '#56B4E9' # Lighter blue
color_female_dart = '#E69F00' # Lighter orange

scaffold_length = 122615

# ==============================
# PANEL A: WGS Depth
# ==============================
ax1.plot(df['start'], df['wgs_m_smooth'], color=color_male_wgs, linewidth=2.5, label='Male (WGS)', zorder=4)
ax1.plot(df['start'], df['wgs_f_smooth'], color=color_female_wgs, linewidth=2.5, alpha=0.8, label='Female (WGS)', zorder=3)

ax1.fill_between(df['start'], df['wgs_m_smooth'], color=color_male_wgs, alpha=0.1, zorder=2)
ax1.fill_between(df['start'], df['wgs_f_smooth'], color=color_female_wgs, alpha=0.1, zorder=1)

ax1.axvspan(58813, 61548, color='grey', alpha=0.15, zorder=0)
max_wgs = max(df['wgs_m_smooth'].max(), df['wgs_f_smooth'].max())
ax1.text((58813+61548)/2, 35, 'Region A\n(Tc1)', color='black', ha='center', fontweight='bold', fontsize=11, zorder=5)

ax1.axvspan(104720, 106467, color='grey', alpha=0.15, zorder=0)
ax1.text((104720+106467)/2, 35, 'Region B\n(LINE)', color='black', ha='center', fontweight='bold', fontsize=11, zorder=5)

ax1.set_ylabel('WGS Depth', fontsize=14, color='black', labelpad=15)
ax1.set_title('A. Whole-Genome WGS Depth (Homomorphic Scaffold Backbone)', fontsize=16, fontweight='bold', pad=15)

ax1.set_ylim(0, 70)
ax1.grid(axis='y', linestyle='--', alpha=0.5, color='grey')
for spine in ['top', 'right']:
    ax1.spines[spine].set_visible(False)
ax1.spines['left'].set_color('black')
ax1.spines['bottom'].set_color('black')
ax1.legend(loc='upper left', frameon=True, facecolor='white', edgecolor='black', fontsize=12)

# ==============================
# PANEL B: DArT-seq Depth
# ==============================
ax2.plot(df['start'], df['dart_m_smooth'], color=color_male_dart, linewidth=2.5, label='Male (DArT-seq)', zorder=4)
ax2.plot(df['start'], df['dart_f_smooth'], color=color_female_dart, linewidth=2.5, alpha=0.8, label='Female (DArT-seq)', zorder=3)

ax2.fill_between(df['start'], df['dart_m_smooth'], color=color_male_dart, alpha=0.1, zorder=2)
ax2.fill_between(df['start'], df['dart_f_smooth'], color=color_female_dart, alpha=0.1, zorder=1)

ax2.axvspan(58813, 61548, color='grey', alpha=0.15, zorder=0)
ax2.axvspan(104720, 106467, color='grey', alpha=0.15, zorder=0)

ax2.set_ylabel('DArT-seq Depth', fontsize=14, color='black', labelpad=15)
ax2.set_xlabel('Genomic Position on Scaffold_2453 (kb)', fontsize=14, color='black', labelpad=15)
ax2.set_title('B. Whole-Genome DArT-seq Depth', fontsize=16, fontweight='bold', pad=15)

ax2.set_ylim(0, 70)
ax2.set_xlim(0, scaffold_length)
ax2.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{int(x/1000)}"))

ax2.grid(axis='y', linestyle='--', alpha=0.5, color='grey')
for spine in ['top', 'right']:
    ax2.spines[spine].set_visible(False)
ax2.spines['left'].set_color('black')
ax2.spines['bottom'].set_color('black')
ax2.tick_params(colors='black', labelsize=12)
ax2.legend(loc='upper left', frameon=True, facecolor='white', edgecolor='black', fontsize=12)

plt.tight_layout()
plt.savefig('Figure_S1_Masterpiece.pdf', dpi=400, bbox_inches='tight', facecolor='white', transparent=False)
plt.savefig('Figure_S1_Masterpiece.png', dpi=400, bbox_inches='tight', facecolor='white', transparent=False)
print("Successfully created dual-panel Figure_S1_Masterpiece.pdf")
