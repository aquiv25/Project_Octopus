"""
Octopus bimaculoides — combined Z-DNA figure.
Top row:    A Genome-wide feature map | B Z-DNABERT | C Z-Hunter  (pie charts)
Bottom row: D ideogram (heatmap on chromosomes)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"]  = "Arial"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
from pathlib import Path

DATA_DIR = Path("/Users/dassagaripova/Downloads/project")

# ══════════════════════════════════════════════════════════════
# PIE DATA
# ══════════════════════════════════════════════════════════════
COLORS = {
    "Intergenic":  "#A5A5A5",
    "Introns":     "#4472C4",
    "Promoters":   "#FFC000",
    "Exons":       "#FF0000",
    "UTR":         "#70AD47",
    "CDS":         "#ED7D31",
    "Downstream":  "#7030A0",
}
CAT_ORDER = ["Intergenic","Introns","Promoters","Exons","UTR","CDS","Downstream"]

pie_data = {
    ("B", "Genome", ""): {
        "Intergenic": 54.50, "Introns": 40.63, "Promoters": 2.06,
        "Exons": 2.27, "UTR": 1.12, "CDS": 1.14, "Downstream": 1.76,
    },
    ("C", "Z-DNABERT", "N = 1,588,238"): {
        "Intergenic": 61.26, "Introns": 29.27, "Promoters": 1.97,
        "Exons": 0.47, "UTR": 0.40, "CDS": 0.07, "Downstream": 1.29,
    },
    ("D", "Z-Hunter", "N = 1,781,183"): {
        "Intergenic": 58.84, "Introns": 38.09, "Promoters": 1.93,
        "Exons": 0.84, "UTR": 0.82, "CDS": 0.02, "Downstream": 1.24,
    },
}

# ══════════════════════════════════════════════════════════════
# IDEOGRAM DATA
# ══════════════════════════════════════════════════════════════
CHR_LENGTHS = [
    199866856,192480835,168055261,147789705,127475568,
    117078234,110520290, 97793126, 96880575, 94596256,
     80862590, 68457376, 65055609, 60575670, 58152616,
     57488271, 57420321, 55529144, 55115435, 47890911,
     40058617, 37711260, 37670884, 36152331, 35487451,
     26644192, 22551867, 18633596, 11972102, 10081836,
]
N_CHR     = 30
CHR_NAMES = [f"Chr{i}" for i in range(1, N_CHR+1)]
MAX_LEN   = max(CHR_LENGTHS)
NC_IDS    = ([f"NC_06{n}.1" for n in range(8981,9000)] +
             [f"NC_06{n}.1" for n in range(9000,9011)])
CHR_MAP   = dict(zip(NC_IDS, CHR_NAMES))

def read_bg(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith(("track","browser","#")) or not line.strip(): continue
            p = line.split()
            if len(p) < 4 or p[0] not in CHR_MAP: continue
            try: rows.append((CHR_MAP[p[0]], int(p[1]), int(p[2]), float(p[3])))
            except: continue
    return pd.DataFrame(rows, columns=["chr","start","end","score"])

print("Reading bedgraphs...")
bert  = read_bg(DATA_DIR/"GCF_001194135.2_ASM119413v2_zdna_thr025.bedgraph.gz")
zhunt = read_bg(DATA_DIR/"zhunter_z-dna.bedgraph")
print(f"  Z-DNABERT {len(bert):,}  Z-Hunter {len(zhunt):,}")

WINDOW = 1_000_000
def density(df):
    res = {}
    for name, length in zip(CHR_NAMES, CHR_LENGTHS):
        sub = df[df["chr"]==name]
        brk = np.arange(0, length+WINDOW, WINDOW)
        vals = []
        for s,e in zip(brk[:-1], brk[1:]):
            e = min(e, length)
            ov = (np.minimum(sub["end"],e)-np.maximum(sub["start"],s)).clip(lower=0).sum() if len(sub) else 0
            vals.append(ov/(e-s)*100)
        res[name] = np.array(vals)
    return res

print("Computing densities...")
bd = density(bert);  vb = np.percentile(np.concatenate(list(bd.values())), 99)
zd = density(zhunt); vz = np.percentile(np.concatenate(list(zd.values())), 99)

cmap_b = LinearSegmentedColormap.from_list("b", ["#FFFFFF","#C6DBEF","#6BAED6","#2171B5","#08306B"])
cmap_z = LinearSegmentedColormap.from_list("z", ["#FFFFFF","#C7E9C0","#74C476","#238B45","#00441B"])

def chr_path(xc, yb, h, hw):
    w=hw; yt=yb+h; ch=w*1.4; N=80
    t1=np.linspace(0,np.pi,N); t2=np.linspace(np.pi,2*np.pi,N)
    ax_=np.concatenate([xc+w*np.cos(t1),[xc+w,xc+w],xc+w*np.cos(t2),[xc-w,xc-w]])
    ay =np.concatenate([yt+ch*np.sin(t1),[yt,yb],    yb+ch*np.sin(t2),[yb,  yt]])
    v=list(zip(ax_,ay))+[(ax_[0],ay[0])]
    c=[MPath.MOVETO]+[MPath.LINETO]*(len(v)-2)+[MPath.CLOSEPOLY]
    return MPath(v,c)

# ══════════════════════════════════════════════════════════════
# NON-OVERLAPPING LABEL PLACEMENT  (exact style of reference)
# ══════════════════════════════════════════════════════════════
def place_pie_labels(ax, wedges, labels, pcts, fs=9.0):
    """
    Two-segment callout:  wedge-edge → kink point → text anchor.
    Left/right separated; labels spread vertically without overlap.
    """
    R_EDGE  = 1.04   # radius where leader starts (just outside wedge)
    R_KINK  = 1.22   # radius of kink/elbow point
    X_R     =  1.72  # fixed x for right-side text
    X_L     = -1.72  # fixed x for left-side text
    MIN_GAP = 0.175  # minimum vertical gap between labels
    FS      = fs

    entries = []
    for w, lbl, pct in zip(wedges, labels, pcts):
        ang = np.deg2rad((w.theta2 + w.theta1) / 2)
        entries.append({"lbl": lbl, "pct": pct, "ang": ang,
                        "cos": np.cos(ang), "sin": np.sin(ang)})

    def spread(ys, gap, lo=-1.30, hi=1.30):
        n = len(ys)
        ys = list(ys)
        for _ in range(800):
            moved = False
            for i in range(n - 1):
                d = ys[i] - ys[i+1]
                if d < gap:
                    push = (gap - d) / 2
                    ys[i]   += push
                    ys[i+1] -= push
                    moved = True
            if not moved:
                break
        # clamp whole group into [lo, hi]
        if ys[0] > hi:
            delta = ys[0] - hi
            ys = [y - delta for y in ys]
        if ys[-1] < lo:
            delta = lo - ys[-1]
            ys = [y + delta for y in ys]
        return ys

    right = sorted([e for e in entries if e["cos"] >= 0],
                   key=lambda e: -e["sin"])
    left  = sorted([e for e in entries if e["cos"] <  0],
                   key=lambda e: -e["sin"])

    for group, x_txt, ha in [(right, X_R, "left"), (left, X_L, "right")]:
        if not group:
            continue
        ys_init = [e["sin"] for e in group]
        ys      = spread(ys_init, MIN_GAP)

        for e, y_lbl in zip(group, ys):
            # segment 1: edge → kink
            x0, y0 = R_EDGE * e["cos"], R_EDGE * e["sin"]
            x1, y1 = R_KINK * e["cos"], R_KINK * e["sin"]
            # segment 2: kink → text anchor (horizontal-ish)
            x2 = x_txt * 0.95
            y2 = y_lbl

            ax.plot([x0, x1, x2], [y0, y1, y2],
                    color="#555", lw=0.65,
                    solid_capstyle="round", solid_joinstyle="round",
                    clip_on=False, zorder=5)

            txt = f"{e['lbl']}, {e['pct']:.2f}%"
            ax.text(x_txt, y_lbl, txt,
                    ha=ha, va="center",
                    fontsize=FS, color="#111",
                    clip_on=False)

# ══════════════════════════════════════════════════════════════
# FIGURE
# ══════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(22, 13), facecolor="white")

outer = gridspec.GridSpec(3, 1, figure=fig,
                          height_ratios=[1.25, 0.06, 2.1],
                          hspace=0.06)

# ── Pie row ───────────────────────────────────────────────────
pie_gs = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[0], wspace=0.30)

for pi, ((letter, title, subtitle), data) in enumerate(pie_data.items()):
    ax = fig.add_subplot(pie_gs[pi])

    sizes  = [data[k] for k in CAT_ORDER]
    total  = sum(sizes)
    pcts   = [s/total*100 for s in sizes]
    cols   = [COLORS[k] for k in CAT_ORDER]

    wedges, _ = ax.pie(pcts, colors=cols, startangle=90,
                       counterclock=False,
                       wedgeprops=dict(linewidth=0.5, edgecolor="white"),
                       radius=1.0)

    place_pie_labels(ax, wedges, CAT_ORDER, pcts, fs=9.0)

    ax.set_xlim(-2.75, 2.75)
    ax.set_ylim(-1.55, 1.80)

    # Title: bold name on line 1, bold N= on line 2
    if subtitle:
        ax.text(0, 1.78, title,    transform=ax.transData,
                ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.text(0, 1.60, subtitle, transform=ax.transData,
                ha="center", va="bottom", fontsize=11, fontweight="bold")
    else:
        ax.text(0, 1.78, title, transform=ax.transData,
                ha="center", va="bottom", fontsize=11, fontweight="bold")

    # Panel letter top-left
    ax.text(-0.04, 1.10, letter, transform=ax.transAxes,
            fontsize=15, fontweight="bold", va="top", ha="right")

# ── Shared legend row ─────────────────────────────────────────
ax_leg = fig.add_subplot(outer[1]); ax_leg.axis("off")
handles = [mpatches.Patch(fc=COLORS[k], ec="white", lw=0.4, label=k)
           for k in CAT_ORDER]
ax_leg.legend(handles=handles, loc="center", ncol=7, fontsize=9.5,
              frameon=False, handlelength=1.3, columnspacing=1.3)

# ── Ideogram ──────────────────────────────────────────────────
bot_gs  = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[2])
ax_ideo = fig.add_subplot(bot_gs[0])
ax_ideo.set_facecolor("white"); ax_ideo.axis("off")
ax_ideo.text(-0.008, 1.005, "A", transform=ax_ideo.transAxes,
             fontsize=15, fontweight="bold", va="bottom", ha="right")

COL_SEP=1.0; HW=0.30
ax_ideo.set_xlim(-1.0, N_CHR*COL_SEP+0.5)
ax_ideo.set_ylim(-0.09*MAX_LEN, MAX_LEN*1.08)

print("Drawing ideogram...")
for idx,(name,length) in enumerate(zip(CHR_NAMES,CHR_LENGTHS)):
    xc=idx*COL_SEP+COL_SEP/2; h=length
    nw=len(bd[name]); wh=h/nw
    fp=chr_path(xc,0,h,HW)
    for wi,val in enumerate(bd[name]):
        yw=wi*wh
        r=mpatches.Rectangle((xc-HW-0.02,yw),HW+0.02,wh+5e4,
                              linewidth=0,facecolor=cmap_b(min(val/vb,1)),zorder=2)
        cp=PathPatch(fp,transform=ax_ideo.transData)
        ax_ideo.add_patch(cp); r.set_clip_path(cp); ax_ideo.add_patch(r)
    for wi,val in enumerate(zd[name]):
        yw=wi*wh
        r=mpatches.Rectangle((xc,yw),HW+0.02,wh+5e4,
                              linewidth=0,facecolor=cmap_z(min(val/vz,1)),zorder=2)
        cp=PathPatch(fp,transform=ax_ideo.transData)
        ax_ideo.add_patch(cp); r.set_clip_path(cp); ax_ideo.add_patch(r)
    ax_ideo.plot([xc,xc],[0,h],color="white",lw=0.35,zorder=9)
    ax_ideo.add_patch(PathPatch(fp,facecolor="none",edgecolor="#222",lw=0.65,zorder=10))
    ax_ideo.text(xc,h+MAX_LEN*0.013,str(idx+1),
                 ha="center",va="bottom",fontsize=6.5,fontweight="bold",color="#111")
    ax_ideo.text(xc,-MAX_LEN*0.011,f"{length/1e6:.0f}",
                 ha="center",va="top",fontsize=4.5,color="#888")

# Scale bar
bx=-0.68
ax_ideo.plot([bx,bx],[0,100e6],color="#333",lw=1.3,zorder=10)
ax_ideo.plot([bx-.04,bx+.04],[0,0],color="#333",lw=0.9,zorder=10)
ax_ideo.plot([bx-.04,bx+.04],[100e6,100e6],color="#333",lw=0.9,zorder=10)
ax_ideo.text(bx-.12,50e6,"100 Mb",ha="right",va="center",
             fontsize=8,color="#333",rotation=90)

# Colorbars
for cmap_i,vm,lbl,yb in [(cmap_b,vb,"Z-DNABERT\n(% bp / 1 Mb)",0.57),
                          (cmap_z,vz,"Z-Hunter\n(% bp / 1 Mb)",  0.32)]:
    ac=fig.add_axes([0.955,yb,0.010,0.19])
    sm=ScalarMappable(cmap=cmap_i); sm.set_array([0,vm])
    cb=fig.colorbar(sm,cax=ac,orientation="vertical")
    cb.set_label(lbl,fontsize=8,labelpad=5)
    cb.outline.set_linewidth(0.4); ac.tick_params(labelsize=7)
    ac.text(2.2,0.0,"Low", transform=ac.transAxes,fontsize=7,va="bottom",color="#555")
    ac.text(2.2,1.0,"High",transform=ac.transAxes,fontsize=7,va="top",   color="#555")

leg2=[mpatches.Patch(fc="#2171B5",label="Z-DNABERT (left)"),
      mpatches.Patch(fc="#238B45",label="Z-Hunter (right)")]
fig.legend(handles=leg2,loc="lower right",bbox_to_anchor=(0.950,0.27),
           fontsize=8.5,frameon=True,framealpha=0.95,edgecolor="#ccc",
           title="Z-DNA predictor",title_fontsize=8.5)

# ══════════════════════════════════════════════════════════════
out_svg = DATA_DIR/"octopus_figure.svg"
out_png = DATA_DIR/"octopus_figure.png"
fig.savefig(out_svg,format="svg",bbox_inches="tight")
fig.savefig(out_png,format="png",bbox_inches="tight",dpi=300)
print(f"\n✅ SVG: {out_svg}\n✅ PNG: {out_png}")
plt.close()
