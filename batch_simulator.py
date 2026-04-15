import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons
from matplotlib.lines import Line2D

# --- 1. Daten generieren (Simulation der Zell-Koordinaten) ---
np.random.seed(42) # Sorgt dafür, dass die "zufälligen" Punkte immer gleich aussehen
n_per_group = 50   # 50 Zellen pro Typ pro Experiment

# Zentren für Rohe Daten: Starker Batch-Effekt (Batch A ist links, Batch B ist rechts)
raw_centers = {
    ('Batch A', 'T-Zellen'): [-5, 3],
    ('Batch A', 'B-Zellen'): [-6, 0],
    ('Batch A', 'Makrophagen'): [-4, -3],
    ('Batch B', 'T-Zellen'): [4, 4],
    ('Batch B', 'B-Zellen'): [5, 1],
    ('Batch B', 'Makrophagen'): [6, -2]
}

# Zentren für korrigierte Daten: Batches überlappen exakt, nur Zelltypen trennen sich
int_centers = {
    ('Batch A', 'T-Zellen'): [0, 5],
    ('Batch A', 'B-Zellen'): [-4, -2],
    ('Batch A', 'Makrophagen'): [4, -2],
    ('Batch B', 'T-Zellen'): [0, 5],
    ('Batch B', 'B-Zellen'): [-4, -2],
    ('Batch B', 'Makrophagen'): [4, -2]
}

data_raw = {}
data_int = {}

# Optik definieren
colors = {'T-Zellen': '#1f77b4', 'B-Zellen': '#ff7f0e', 'Makrophagen': '#2ca02c'}
markers = {'Batch A': 'o', 'Batch B': 'X'}

# Punktewolken berechnen (Zentrum + zufälliges Rauschen)
for batch in ['Batch A', 'Batch B']:
    for ctype in ['T-Zellen', 'B-Zellen', 'Makrophagen']:
        data_raw[(batch, ctype)] = raw_centers[(batch, ctype)] + np.random.randn(n_per_group, 2) * 1.2
        data_int[(batch, ctype)] = int_centers[(batch, ctype)] + np.random.randn(n_per_group, 2) * 1.2

# --- 2. Plot aufbauen ---
fig, ax = plt.subplots(figsize=(11, 7)) # Breite leicht erhöht für mehr Platz
fig.canvas.manager.set_window_title('Batch Effect Integration Simulator')

# Wir zwingen den Hauptgraphen weiter nach rechts (Start bei 35% statt 30%)
plt.subplots_adjust(left=0.35, right=0.95, top=0.9, bottom=0.1) 

# Alle Punktewolken initial in den Plot laden (Zustand: Rohdaten)
scatters = {}
for batch in ['Batch A', 'Batch B']:
    for ctype in ['T-Zellen', 'B-Zellen', 'Makrophagen']:
        pts = data_raw[(batch, ctype)]
        scatters[(batch, ctype)] = ax.scatter(
            pts[:, 0], pts[:, 1], 
            c=colors[ctype], 
            marker=markers[batch], 
            alpha=0.7, 
            edgecolor='white',
            s=80 
        )

# Achsen formatieren
ax.set_title("Single-Cell Simulation: Batch-Effekt vs. Integration", fontweight='bold')
ax.set_xlabel("UMAP / PCA Dimension 1")
ax.set_ylabel("UMAP / PCA Dimension 2")
ax.set_xlim(-10, 10)
ax.set_ylim(-8, 8)
ax.grid(True, linestyle='--', alpha=0.3)

# Eine schöne, aufgeräumte Legende bauen
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Batch A (Exp. 1)', markerfacecolor='gray', markersize=10),
    Line2D([0], [0], marker='X', color='w', label='Batch B (Exp. 2)', markerfacecolor='gray', markersize=10),
    Line2D([0], [0], color='w', label=''), # Abstandshalter
    Line2D([0], [0], marker='o', color='w', label='T-Zellen', markerfacecolor=colors['T-Zellen'], markersize=10),
    Line2D([0], [0], marker='o', color='w', label='B-Zellen', markerfacecolor=colors['B-Zellen'], markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Makrophagen', markerfacecolor=colors['Makrophagen'], markersize=10)
]

# Legende GANZ OBEN LINKS aufhängen (y=0.95), sodass sie sicher nach unten wächst
fig.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0.02, 0.95), frameon=True, framealpha=0.9)

# --- 3. Interaktivität (Radio Buttons) ---
# Box für die Buttons etwas weiter unten platzieren (y=0.50), damit Abstand zur Legende ist
ax_radio = plt.axes([0.02, 0.50, 0.25, 0.15])
ax_radio.set_facecolor('#f0f0f0')
radio = RadioButtons(ax_radio, ('Rohe Daten\n(Batch-Effekt)', 'Korrigierte Daten\n(Integriert)'))

# Funktion, die ausgeführt wird, wenn man einen Button klickt
def update(label):
    for batch in ['Batch A', 'Batch B']:
        for ctype in ['T-Zellen', 'B-Zellen', 'Makrophagen']:
            if label.startswith('Rohe'):
                # Setze die Punkte auf die unkorrigierten Koordinaten zurück
                scatters[(batch, ctype)].set_offsets(data_raw[(batch, ctype)])
            else:
                # Setze die Punkte auf die korrigierten, überlappenden Koordinaten
                scatters[(batch, ctype)].set_offsets(data_int[(batch, ctype)])
    # Bild neu zeichnen
    fig.canvas.draw_idle()

# Die Funktion an die Buttons koppeln
radio.on_clicked(update)

plt.show()