import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# --- 1. Daten und Mathematik ---
# X-Achse (Genomische Position 0 bis 1000)
x = np.linspace(0, 1000, 1000)

# Positionen der offenen Regionen
pos_enh_A = 250
pos_enh_B = 500
pos_gen_X = 800
peak_width = 15

# Funktion für die Form der Peaks (Gauß-Glocke)
def gauss(x_vals, center, height, width):
    return height * np.exp(-((x_vals - center)**2) / (2 * width**2))

# Die Profile der drei Einzelzellen generieren
prof_yellow = gauss(x, pos_enh_A, 1.0, peak_width) + gauss(x, pos_gen_X, 1.0, peak_width)
prof_green = gauss(x, pos_enh_B, 1.0, peak_width) + gauss(x, pos_gen_X, 1.0, peak_width)
prof_blue = gauss(x, pos_enh_A, 0.3, peak_width) + gauss(x, pos_gen_X, 1.0, peak_width)


# --- 2. Das Plot-Layout (Matplotlib) ---
# Erstelle ein Fenster mit 4 untereinanderliegenden Plots (Subplots)
fig, (ax_bulk, ax_y, ax_g, ax_b) = plt.subplots(4, 1, figsize=(10, 8), sharex=True)
plt.subplots_adjust(bottom=0.3, hspace=0.3) # Platz unten für Slider machen
fig.canvas.manager.set_window_title('Bulk vs scATAC-seq Simulator')

# Farben definieren
color_bulk = 'gray'
color_y = '#fcc419'
color_g = '#40c057'
color_b = '#339af0'

# Die statischen Einzelzell-Spuren zeichnen (als gefüllte Flächen)
ax_y.fill_between(x, 0, prof_yellow, color=color_y, alpha=0.8)
ax_g.fill_between(x, 0, prof_green, color=color_g, alpha=0.8)
ax_b.fill_between(x, 0, prof_blue, color=color_b, alpha=0.8)

# Achsenbeschriftungen und Formatierung
axes = [ax_bulk, ax_y, ax_g, ax_b]
titles = ['Bulk Gewebe (Mischung)', 'Zelltyp Gelb', 'Zelltyp Grün', 'Zelltyp Blau']

for ax, title in zip(axes, titles):
    ax.set_ylabel(title, rotation=0, labelpad=60, ha='center', va='center', fontweight='bold')
    ax.set_yticks([]) # Y-Zahlen ausblenden
    ax.set_ylim(0, 1.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Vertikale gestrichelte Linien für die Enhancer/Promotor Zonen
    for pos, label in zip([pos_enh_A, pos_enh_B, pos_gen_X], ['Enhancer A', 'Enhancer B', 'Gen X Promotor']):
        ax.axvline(pos, color='black', linestyle='--', alpha=0.2)
        if ax == ax_b: # Nur beim untersten Plot die Namen hinschreiben
            ax.text(pos, -0.4, label, ha='center', va='top', fontsize=10, color='gray')

ax_b.set_xticks([]) # X-Zahlen ausblenden

# --- 3. Die Interaktivität (Slider & Update-Funktion) ---
# Achsenbereiche für die Slider definieren (Position: Links, Unten, Breite, Höhe)
ax_slider_y = plt.axes([0.2, 0.15, 0.6, 0.03])
ax_slider_g = plt.axes([0.2, 0.10, 0.6, 0.03])
ax_slider_b = plt.axes([0.2, 0.05, 0.6, 0.03])

# Slider erstellen
slider_y = Slider(ax_slider_y, 'Anzahl Gelb', 0, 100, valinit=33, valstep=1, color=color_y)
slider_g = Slider(ax_slider_g, 'Anzahl Grün', 0, 100, valinit=33, valstep=1, color=color_g)
slider_b = Slider(ax_slider_b, 'Anzahl Blau', 0, 100, valinit=33, valstep=1, color=color_b)

# Diese Funktion wird aufgerufen, wenn ein Slider bewegt wird
def update(val):
    w_y = slider_y.val
    w_g = slider_g.val
    w_b = slider_b.val
    total = w_y + w_g + w_b
    
    if total > 0:
        # Die gewichtete Summe berechnen
        prof_bulk = (prof_yellow * w_y + prof_green * w_g + prof_blue * w_b) / total
    else:
        prof_bulk = np.zeros_like(x)
    
    # Alten Bulk-Plot löschen und neuen zeichnen
    for coll in ax_bulk.collections:
        coll.remove()
    ax_bulk.fill_between(x, 0, prof_bulk, color=color_bulk, alpha=0.8)
    
    # Y-Achse dynamisch anpassen, damit der höchste Peak gut sichtbar bleibt
    max_val = max(prof_bulk.max() * 1.1, 0.1)
    ax_bulk.set_ylim(0, max_val)
    
    # Das Bild aktualisieren
    fig.canvas.draw_idle()

# Die Update-Funktion an die Slider koppeln
slider_y.on_changed(update)
slider_g.on_changed(update)
slider_b.on_changed(update)

# Die Funktion einmal manuell aufrufen, um das Startbild zu generieren
update(0)

# Den Plot anzeigen
plt.show()