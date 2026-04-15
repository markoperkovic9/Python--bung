import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# --- 1. Daten generieren (Die 3D Swiss Roll) ---
np.random.seed(42)
n_samples = 1500
noise = 0.05

# Parametrische Gleichungen für die Swiss Roll
# t bestimmt den Winkel und den Radius der Rolle (wie weit wir aufgerollt sind)
t = 1.5 * np.pi * (1 + 2 * np.random.rand(n_samples))
x_rolled = t * np.cos(t)
y_rolled = 21 * np.random.rand(n_samples) # Die Höhe der Rolle
z_rolled = t * np.sin(t)

# Die 3D Koordinaten zusammenfügen und etwas Rauschen hinzufügen
X_rolled = np.vstack((x_rolled, y_rolled, z_rolled)).T + noise * np.random.randn(n_samples, 3)

# --- 2. Definition des entrollten 2D-Zustands ---
# Die wahre entrollte Struktur (das Ziel von UMAP)
# Wir mappen 't' auf die horizontale Achse (die abgerollte Länge)
# und behalten 'y' als Höhe. Die Tiefe 'z' wird flach (0).
x_unrolled = t
y_unrolled = y_rolled
z_unrolled = np.zeros(n_samples)

X_unrolled = np.vstack((x_unrolled, y_unrolled, z_unrolled)).T

# Wir zentrieren den entrollten Zustand für eine bessere Visualisierung
X_unrolled[:, 0] -= X_unrolled[:, 0].mean()


# --- 3. Färbung definieren (Regenbogen entlang der Rolle) ---
# Das ist entscheidend: Wir färben basierend auf dem Wert von 't'.
# So sehen wir, ob die Farbreihenfolge beim Entrollen erhalten bleibt.
colors = plt.cm.jet((t - t.min()) / (t.max() - t.min()))


# --- 4. Plot aufbauen (Matplotlib) ---
fig = plt.figure(figsize=(10, 8))
fig.canvas.manager.set_window_title('Single-Cell UMAP Manifold Unrolling Simulator')

# Standard 3D Achsensystem erstellen
ax = fig.add_subplot(111, projection='3d')
plt.subplots_adjust(bottom=0.25) # Platz unten für Slider machen

# Initialer Plot: 3D Aufgerollt
scatter = ax.scatter(X_rolled[:,0], X_rolled[:,1], X_rolled[:,2], c=colors, s=15, alpha=0.7, edgecolors='none')

# Statische Achsengrenzen setzen (wichtig für fließende Animation)
all_X = np.concatenate((X_rolled, X_unrolled))
ax.set_xlim(all_X[:,0].min() - 1, all_X[:,0].max() + 1)
ax.set_ylim(all_X[:,1].min() - 1, all_X[:,1].max() + 1)
ax.set_zlim(all_X[:,2].min() - 1, all_X[:,2].max() + 1)

# Achsen beschriften und formatieren
ax.set_title("UMAP Simulation: Manifold Unrolling\nBeobachte, wie die lokale Struktur erhalten bleibt!", fontweight='bold', pad=20)
ax.set_xlabel("Gen-Expression 1")
ax.set_ylabel("Gen-Expression 2")
ax.set_zlabel("Gen-Expression 3")

# Visueller Hinweis zur Drehbarkeit
fig.text(0.5, 0.95, "Drehe die Rolle mit der Maus!", ha='center', fontweight='bold', fontsize=12, color='white')


# --- 5. Interaktivität (Slider & Animation) ---
# Achsenbereich für den Slider definieren
ax_slider = plt.axes([0.2, 0.1, 0.6, 0.03])
# Slider erstellen (Wert 0.0 = Aufgerollt, 1.0 = Entrollt)
# Farbwahl für den Slider angepasst an meine Simulation
slider = Slider(ax_slider, 'UMAP-Entrollung (3D → 2D)', 0.0, 1.0, valinit=0.0, color='royalblue')

# Diese Funktion wird ausgeführt, wenn der Slider bewegt wird
def update(val):
    alpha = slider.val
    # Die lineare Interpolation: Wir berechnen die Zwischenstationen
    # zwischen aufgerollt (alpha=0) und entrollt (alpha=1)
    X_curr = (1 - alpha) * X_rolled + alpha * X_unrolled
    
    # Die Koordinaten im 3D-Plot aktualisieren (Schattenwurf / Projektion)
    scatter._offsets3d = (X_curr[:, 0], X_curr[:, 1], X_curr[:, 2])
    fig.canvas.draw_idle()

# Die Update-Funktion an den Slider koppeln
slider.on_changed(update)

# Den Plot anzeigen
plt.show()