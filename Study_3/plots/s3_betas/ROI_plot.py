import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from nilearn import plotting
import io
from PIL import Image

# ── Coordinates ─────────────────────────────────────────────
coords = [
    [-44, 22, -8],   # lIFG
    [56, 24, 18],    # rIFG
    [-62, -34, -6],  # lMTG
    [62, -30, -6],   # rMTG
    [-44, -56, 40],  # lTPJ
    [50, -56, 30],   # rTPJ
    [-4, 42, 28],    # lmPFC
    [6, 42, 36],     # rmPFC
    [-4, -66, 32],   # lPCun
]

region_groups = {
    "IFG"  : {"indices": [0, 1], "color": "#0072B2"},
    "MTG"  : {"indices": [2, 3], "color": "#E69F00"},
    "TPJ"  : {"indices": [4, 5], "color": "#009E73"},
    "mPFC" : {"indices": [6, 7], "color": "#CC79A7"},
    "PCun" : {"indices": [8],    "color": "#D55E00"},
}

legend_patches = [
    mpatches.Patch(color=info["color"], label=region)
    for region, info in region_groups.items()
]

# ── Shared font sizes (mirror the R theme) ───────────────────
LABEL_SIZE  = 9    # strip labels / view labels  → strip.text = 9
LEGEND_SIZE = 9    # legend text / title          → legend.text/title = 9
TITLE_SIZE  = 10   # panel title                  → plot.title = 10

# ── Render each projection into a buffer ────────────────────
def render_view(mode):
    fig_tmp = plt.figure(figsize=(3, 2.5), facecolor="white")
    display = plotting.plot_glass_brain(
        None,
        display_mode=mode,
        colorbar=False,
        plot_abs=False,
        figure=fig_tmp,
    )
    for region, info in region_groups.items():
        region_coords = [coords[i] for i in info["indices"]]
        display.add_markers(
            region_coords,
            marker_color=info["color"],
            marker_size=90,
        )
    buf = io.BytesIO()
    fig_tmp.savefig(buf, format="png", dpi=150,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig_tmp)
    buf.seek(0)
    return Image.open(buf).copy()

view_modes  = ["l", "y", "r", "z"]
view_labels = ["Sagittal L", "Coronal", "Sagittal R", "Axial"]
images      = [render_view(m) for m in view_modes]

# ── Stack 4 views vertically ────────────────────────────────
fig, axes = plt.subplots(4, 1, figsize=(3, 11), facecolor="white")
fig.subplots_adjust(hspace=0.05, top=0.93)

# # Panel title 
# fig.text(
#     0.5, 0.965,
#     "ROI anatomical localisation",
#     ha="center", va="top",
#     fontsize=TITLE_SIZE,
#     fontweight="normal",
#     color="black"
# )

for ax, img, label in zip(axes, images, view_labels):
    ax.imshow(img)
    ax.axis("off")
    # View label — matches strip.text size = 9, bold

# ── Legend at bottom — matches legend.text/title size = 9 ───
fig.legend(
    handles=legend_patches,
    loc="lower center",
    ncol=2,
    frameon=False,
    fontsize=LEGEND_SIZE,
    bbox_to_anchor=(0.5, 0.0),
    title="ROIs",
    title_fontsize=LEGEND_SIZE,
)

plt.savefig(
    "roi_glass_brain_vertical.png",
    dpi=600,
    bbox_inches="tight",
    facecolor="white"
)
plt.show()