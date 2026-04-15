import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# --- 1. Daten generieren (4 Cluster in 3D) ---
np.random.seed(42)
n_cells = 40

# Wir erschaffen 4 fiktive Zelltypen in einem 3D-Gen-Raum (Gen 1, Gen 2, Gen 3)
# Die Wolke ist absichtlich länglich / schräg im Raum verteilt
c1 = np.random.randn(n_cells, 3) * 0.8 + [1, 1, 1]
c2 = np.random.randn(n_cells, 3) * 0.8 + [4, 3, 2]
c3 = np.random.randn(n_cells, 3) * 0.8 + [7, 6, 4]
c4 = np.random.randn(n_cells, 3) * 0.8 + [3, 7, 5]

# Alle Daten zusammenfügen
X_original = np.vstack([c1, c2, c3, c4])

# Farben für die Zelltypen
colors = ['#ff7f0e']*n_cells + ['#2ca02c']*n_cells + ['#1f77b4']*n_cells + ['#e377c2']*n_cells


# --- 2. Die PCA Mathematik (ohne externe Bibliotheken, nur NumPy) ---
# Schritt 1: Daten zentrieren (Mittelpunkt auf 0,0,0 schieben)
X_mean = np.mean(X_original, axis=0)
X_centered = X_original - X_mean

# Schritt 2: Hauptkomponenten (PCs) berechnen mittels Singulärwertzerlegung (SVD)
U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
PC = Vt[:2, :] # Wir nehmen nur die zwei wichtigsten Komponenten (PC1 und PC2)

# Schritt 3: Die 3D-Punkte auf die 2D-Ebene projizieren (Schattenwurf)
Z_2d = X_centered @ PC.T 

# Schritt 4: Die 2D-Punkte zurück in den 3D-Raum übersetzen (damit sie flach auf der Ebene liegen)
X_projected = (Z_2d @ PC) + X_mean


# --- 3. Den Plot aufbauen ---
fig = plt.figure(figsize=(10, 8))
fig.canvas.manager.set_window_title('PCA 3D Projection Simulator')

# 3D Achsensystem erstellen
ax = fig.add_subplot(111, projection='3d')
plt.subplots_adjust(bottom=0.25) # Platz für Slider machen

# Die Punkte zeichnen (Startzustand: 3D)
scat = ax.scatter(X_original[:,0], X_original[:,1], X_original[:,2], 
                  c=colors, s=50, alpha=0.8, edgecolor='w', linewidth=0.5)

# --- 4. Die PC-Ebene visualisieren (als halbtransparentes Gitter) ---
# Ein Gitternetz erstellen, das der PC1/PC2 Ebene entspricht
grid_range = np.linspace(-6, 6, 10)
grid_x, grid_y = np.meshgrid(grid_range, grid_range)
plane_3d = np.zeros((10, 10, 3))

# Gitter-Punkte in den 3D-Raum drehen
for i in range(10):
    for j in range(10):
        plane_3d[i, j, :] = grid_x[i, j] * PC[0] + grid_y[i, j] * PC[1] + X_mean

# Die Ebene in den Raum zeichnen
ax.plot_surface(plane_3d[:,:,0], plane_3d[:,:,1], plane_3d[:,:,2], 
                alpha=0.15, color='gray', edgecolor='black', linewidth=0.5)


# Achsen beschriften und formatieren
ax.set_title("Dimensionalitätsreduktion: PCA", fontweight='bold', pad=20)
ax.set_xlabel("Gen 1 Expression")
ax.set_ylabel("Gen 2 Expression")
ax.set_zlabel("Gen 3 Expression")

# --- 5. Interaktivität (Slider & Animation) ---
ax_slider = plt.axes([0.2, 0.1, 0.6, 0.03])
slider = Slider(ax_slider, 'Projektion\n(3D -> 2D)', 0.0, 1.0, valinit=0.0, color='royalblue')

def update(val):
    t = slider.val
    # Die Interpolation: Wir berechnen die Zwischenstationen zwischen 3D und "plattgedrückt"
    X_current = (1 - t) * X_original + t * X_projected
    
    # Koordinaten im 3D-Plot aktualisieren
    scat._offsets3d = (X_current[:,0], X_current[:,1], X_current[:,2])
    fig.canvas.draw_idle()

slider.on_changed(update)

plt.show()