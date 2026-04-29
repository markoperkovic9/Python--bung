import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Set up the figure and axes for a 1x4 panel plot
fig, axes = plt.subplots(1, 4, figsize=(20, 5))
fig.suptitle("Step-by-Step Mechanism: Enhancer-Promoter Looping", fontsize=20, fontweight='bold', y=1.05)

# Shared settings for all panels
for ax in axes:
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off') # Hide axes for a cleaner diagram

# ---------------------------------------------------------
# Step 1: Linear DNA & Regulatory Elements
# ---------------------------------------------------------
ax1 = axes[0]
ax1.set_title("1. Linear DNA", fontsize=14, fontweight='bold')

# DNA strand (straight line)
ax1.plot([1, 9], [5, 5], color='black', linewidth=3, zorder=1)

# Enhancer (Left)
enhancer = patches.Rectangle((2, 4.8), 1.5, 0.4, color='orange', zorder=2)
ax1.add_patch(enhancer)
ax1.text(2.75, 4.3, "Enhancer", ha='center', fontsize=12, color='orange', fontweight='bold')

# Promoter (Right)
promoter = patches.Rectangle((7, 4.8), 1, 0.4, color='green', zorder=2)
ax1.add_patch(promoter)
ax1.text(7.5, 4.3, "Promoter\n(TNNT2)", ha='center', fontsize=12, color='green', fontweight='bold')

# ---------------------------------------------------------
# Step 2: Transcription Factor Binding
# ---------------------------------------------------------
ax2 = axes[1]
ax2.set_title("2. TF & Co-activator Binding", fontsize=14, fontweight='bold')

# DNA and Elements
ax2.plot([1, 9], [5, 5], color='black', linewidth=3, zorder=1)
ax2.add_patch(patches.Rectangle((2, 4.8), 1.5, 0.4, color='orange', zorder=2))
ax2.add_patch(patches.Rectangle((7, 4.8), 1, 0.4, color='green', zorder=2))

# Transcription Factor (Circle) binding to Enhancer
tf = patches.Circle((2.75, 5.5), 0.4, color='dodgerblue', zorder=3)
ax2.add_patch(tf)
ax2.text(2.75, 5.5, "TF", ha='center', va='center', fontsize=10, color='white', fontweight='bold')

# ---------------------------------------------------------
# Step 3: Chromatin Looping Initiated
# ---------------------------------------------------------
ax3 = axes[2]
ax3.set_title("3. 3D Chromatin Looping", fontsize=14, fontweight='bold')

# Bending DNA strand (using a simple parabola)
x = np.linspace(2.75, 7.5, 100)
y = 0.4 * (x - 5.125)**2 + 3 # Parabola opening upwards
ax3.plot(x, y, color='black', linewidth=3, zorder=1)
ax3.plot([1, 2.75], [5.25, y[0]], color='black', linewidth=3, zorder=1) # Left tail
ax3.plot([7.5, 9], [y[-1], 5.25], color='black', linewidth=3, zorder=1) # Right tail

# Elements mapped to the bend
ax3.add_patch(patches.Rectangle((2.3, 5), 0.8, 0.4, angle=-30, color='orange', zorder=2))
ax3.add_patch(patches.Rectangle((7.1, 5), 0.8, 0.4, angle=30, color='green', zorder=2))

# TF
tf_loop = patches.Circle((2.6, 5.7), 0.4, color='dodgerblue', zorder=3)
ax3.add_patch(tf_loop)

# Cohesin ring (extruding the loop)
cohesin = patches.Ellipse((5.125, 3), 1.5, 0.8, fill=False, edgecolor='purple', linewidth=3, zorder=4)
ax3.add_patch(cohesin)
ax3.text(5.125, 2.2, "Cohesin", ha='center', fontsize=10, color='purple', fontweight='bold')

# ---------------------------------------------------------
# Step 4: Active Transcription via Mediator Complex
# ---------------------------------------------------------
ax4 = axes[3]
ax4.set_title("4. Active Transcription", fontsize=14, fontweight='bold')

# Tight loop (Enhancer and Promoter touching)
theta = np.linspace(0, np.pi, 100)
x_loop = 5 + 1.5 * np.cos(theta)
y_loop = 4 - 2 * np.sin(theta)
ax4.plot(x_loop, y_loop, color='black', linewidth=3, zorder=1)
ax4.plot([1, 3.5], [4, 4], color='black', linewidth=3, zorder=1) # Left tail
ax4.plot([6.5, 9], [4, 4], color='black', linewidth=3, zorder=1) # Right tail

# Enhancer and Promoter adjacent
ax4.add_patch(patches.Rectangle((3.5, 3.8), 0.8, 0.4, color='orange', zorder=2))
ax4.add_patch(patches.Rectangle((5.7, 3.8), 0.8, 0.4, color='green', zorder=2))

# TF
ax4.add_patch(patches.Circle((3.9, 4.5), 0.4, color='dodgerblue', zorder=3))

# RNA Polymerase II
pol2 = patches.Ellipse((6.1, 4.5), 1.2, 0.8, color='firebrick', zorder=3)
ax4.add_patch(pol2)
ax4.text(6.1, 4.5, "Pol II", ha='center', va='center', fontsize=10, color='white', fontweight='bold')

# Mediator Complex (Bridging TF and Pol II)
mediator = patches.Ellipse((5, 4.8), 1.8, 1, color='gold', zorder=2, alpha=0.8)
ax4.add_patch(mediator)
ax4.text(5, 4.8, "Mediator", ha='center', va='center', fontsize=10, color='black', fontweight='bold')

# mRNA transcript (squiggly line coming off Pol II)
x_rna = np.linspace(6.5, 8.5, 50)
y_rna = 5.5 + 0.2 * np.sin(10 * x_rna) + 0.5 * (x_rna - 6.5)
ax4.plot(x_rna, y_rna, color='red', linewidth=2, zorder=1)
ax4.text(8.7, 6.7, "mRNA", ha='left', fontsize=12, color='red', fontweight='bold')

plt.tight_layout()
plt.show()