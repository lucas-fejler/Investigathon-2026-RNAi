"""Plot top-decile Jaccard agreement from the benchmark comparison table."""

import base64
import csv
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
NAMES = {
    "DM_ATP5F1B": "ATP sintasa β", "DM_BURS_A": "Burs-α", "DM_ETH": "ETH",
    "DM_RPN6": "Rpn6", "DM_VATP_B": "V-ATPasa B", "EPILAB_ATPaseD_Dm": "V-ATPasa D",
    "EPILAB_BicC_Dm": "BicC", "EPILAB_CP_MRFV": "CP (MRFV)",
    "EPILAB_MP_MRFV": "MP putativa (MRFV)",
}


def main():
    """Generate an SVG with a fixed zero-to-one scale and masked missing values.

    Args:
        None.

    Returns:
        None.
    """
    source = ROOT / "benchmark/resultados/comparacion_tres.tsv"
    with source.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    keys = ("jaccard_sifi_dsrna", "jaccard_sifi_siren", "jaccard_dsrna_siren")
    matrix = np.array([
        [float(row[key]) if row[key] not in ("sin datos", "NA") else np.nan
         for key in keys] for row in rows
    ])
    matplotlib.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                               "svg.hashsalt": "investigathon-rnai-top10"})
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("#edf0f0")
    fig, ax = plt.subplots(figsize=(8.4, 4.8), layout="constrained")
    image = ax.imshow(np.ma.masked_invalid(matrix), vmin=0, vmax=1,
                      cmap=cmap, aspect="auto")
    ax.set_yticks(range(len(rows)), [NAMES[row["gen"]] for row in rows])
    ax.set_xticks(range(3), ["si-Fi21 /\ndsRNAEngineer", "si-Fi21 /\nSIREN",
                           "dsRNAEngineer /\nSIREN"])
    ax.tick_params(axis="both", length=0, pad=9)
    ax.set_xticks(np.arange(-0.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    for (row, column), value in np.ndenumerate(matrix):
        label = "N/E" if np.isnan(value) else f"{value:.2f}".replace(".", ",")
        ax.text(column, row, label, ha="center", va="center",
                color="#142d33", fontsize=10)
    for spine in ax.spines.values():
        spine.set_visible(False)
    bar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.035,
                       ticks=[0, 0.25, 0.5, 0.75, 1])
    bar.set_label("Índice de Jaccard")
    bar.ax.yaxis.set_major_formatter(
        FuncFormatter(lambda value, position: f"{value:.2f}".replace(".", ",")))
    bar.outline.set_visible(False)
    destination = Path(__file__).resolve().parent / "figuras/concordancia_top10.svg"
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, metadata={"Date": None}, bbox_inches="tight")
    plt.close(fig)
    svg = destination.read_text(encoding="utf-8")
    destination.write_text(
        "\n".join(line.rstrip() for line in svg.splitlines()) + "\n",
        encoding="utf-8", newline="\n")
    report = destination.parent.parent / "benchmark_rnai.html"
    if report.is_file():
        encoded = base64.b64encode(destination.read_bytes()).decode("ascii")
        uri = "data:image/svg+xml;base64," + encoded
        document, replacements = re.subn(
            r'(<img\s+id="concordancia-top10"\s+src=")[^"]*(")',
            lambda match: match.group(1) + uri + match.group(2),
            report.read_text(encoding="utf-8"), count=1)
        if replacements != 1:
            raise ValueError("The report must contain its concordance figure")
        report.write_text(document, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
