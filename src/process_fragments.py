from pathlib import Path
from collections import defaultdict
import time
import gc

import numpy as np
import pandas as pd


# ============================================================
# PROJECT SETTINGS
# ============================================================

BASE_DIR = Path(__file__).resolve().parent

CANCER_DIR = BASE_DIR / "breast_cancer"
HEALTHY_DIR = BASE_DIR / "healthy"
OUTPUT_DIR = BASE_DIR / "processed"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

COLUMN_NAMES = ["chrom", "start", "stop", "mapq"]

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}

CHUNK_SIZE = 1_000_000
BIN_SIZE = 5_000_000

SHORT_MIN = 100
SHORT_MAX = 150

LONG_MIN = 151
LONG_MAX = 220

MAPQ_MIN = 30
MAX_FRAGMENT_LENGTH = 1000


# ============================================================
# PROCESS ONE SAMPLE
# ============================================================

def process_sample(file_path, group):

    sample_id = file_path.name.split(".")[0]

    qc_output = OUTPUT_DIR / f"{sample_id}_qc_summary.csv"
    global_output = OUTPUT_DIR / f"{sample_id}_global_length_counts.csv"
    regional_output = OUTPUT_DIR / f"{sample_id}_regional_counts.csv"

    # Resume protection
    if (
        qc_output.exists()
        and global_output.exists()
        and regional_output.exists()
    ):
        print(f"{sample_id}: already completed — skipping")
        return

    print("\n" + "=" * 70)
    print(f"Processing: {sample_id}")
    print(f"Group: {group}")
    print(f"File: {file_path.name}")
    print("=" * 70)

    start_time = time.time()

    length_counts = np.zeros(
        MAX_FRAGMENT_LENGTH + 1,
        dtype=np.int64
    )

    regional_counts = defaultdict(
        lambda: {"short": 0, "long": 0}
    )

    total_fragments = 0
    mapq30_fragments = 0
    usable_autosomal_fragments = 0
    regional_fragments = 0
    chunk_number = 0

    reader = pd.read_csv(
        file_path,
        sep="\t",
        header=None,
        names=COLUMN_NAMES,
        usecols=[0, 1, 2, 3],
        compression="gzip",
        chunksize=CHUNK_SIZE,
        dtype={
            "chrom": "string",
            "start": "int32",
            "stop": "int32",
            "mapq": "int16",
        },
    )

    for chunk in reader:

        chunk_number += 1
        total_fragments += len(chunk)

        chrom = chunk["chrom"]

        starts = chunk["start"].to_numpy()
        stops = chunk["stop"].to_numpy()
        mapq = chunk["mapq"].to_numpy()

        lengths = stops - starts

        mapq_mask = mapq >= MAPQ_MIN
        autosome_mask = chrom.isin(AUTOSOMES).to_numpy()

        valid_mask = (
            mapq_mask
            & autosome_mask
            & (lengths > 0)
            & (lengths <= MAX_FRAGMENT_LENGTH)
        )

        mapq30_fragments += int(mapq_mask.sum())

        usable_autosomal_fragments += int(
            valid_mask.sum()
        )

        # ----------------------------------------------------
        # GLOBAL FRAGMENT LENGTH DISTRIBUTION
        # ----------------------------------------------------

        length_counts += np.bincount(
            lengths[valid_mask],
            minlength=MAX_FRAGMENT_LENGTH + 1
        )

        # ----------------------------------------------------
        # REGIONAL SHORT/LONG FRAGMENTS
        # ----------------------------------------------------

        regional_mask = (
            valid_mask
            & (lengths >= SHORT_MIN)
            & (lengths <= LONG_MAX)
        )

        regional_fragments += int(
            regional_mask.sum()
        )

        if regional_mask.any():

            selected_lengths = lengths[regional_mask]

            midpoints = (
                starts[regional_mask].astype(np.int64)
                + stops[regional_mask].astype(np.int64)
            ) // 2

            bin_starts = (
                midpoints // BIN_SIZE
            ) * BIN_SIZE

            size_class = np.where(
                selected_lengths <= SHORT_MAX,
                "short",
                "long"
            )

            regional_chunk = pd.DataFrame({
                "chrom": chrom.loc[
                    regional_mask
                ].to_numpy(),
                "bin_start": bin_starts,
                "size_class": size_class,
            })

            grouped = (
                regional_chunk
                .groupby(
                    ["chrom", "bin_start", "size_class"],
                    observed=True,
                )
                .size()
            )

            for (
                chromosome,
                bin_start,
                fragment_class,
            ), count in grouped.items():

                key = (
                    chromosome,
                    int(bin_start)
                )

                regional_counts[key][
                    fragment_class
                ] += int(count)

        if chunk_number % 10 == 0:

            elapsed = (
                time.time() - start_time
            ) / 60

            print(
                f"  {total_fragments:,} fragments scanned "
                f"| {elapsed:.1f} min"
            )

        del chunk
        gc.collect()

    # ========================================================
    # GLOBAL OUTPUT
    # ========================================================

    length_table = pd.DataFrame({
        "sample_id": sample_id,
        "group": group,
        "fragment_length": np.arange(
            MAX_FRAGMENT_LENGTH + 1
        ),
        "count": length_counts,
    })

    length_table = length_table.loc[
        length_table["count"] > 0
    ].reset_index(drop=True)

    # ========================================================
    # REGIONAL OUTPUT
    # ========================================================

    regional_rows = []

    for (
        chromosome,
        bin_start,
    ), counts in regional_counts.items():

        short_count = counts["short"]
        long_count = counts["long"]

        regional_rows.append({
            "sample_id": sample_id,
            "group": group,
            "chrom": chromosome,
            "bin_start": bin_start,
            "bin_end": bin_start + BIN_SIZE,
            "short_count": short_count,
            "long_count": long_count,
            "total_count":
                short_count + long_count,
            "short_long_ratio":
                short_count / long_count
                if long_count > 0
                else np.nan,
        })

    regional_table = pd.DataFrame(
        regional_rows
    )

    if not regional_table.empty:

        regional_table["chrom_number"] = (
            regional_table["chrom"]
            .str.replace(
                "chr",
                "",
                regex=False
            )
            .astype(int)
        )

        regional_table = (
            regional_table
            .sort_values(
                ["chrom_number", "bin_start"]
            )
            .drop(columns="chrom_number")
            .reset_index(drop=True)
        )

    # ========================================================
    # QC OUTPUT
    # ========================================================

    elapsed_minutes = (
        time.time() - start_time
    ) / 60

    summary_table = pd.DataFrame([{
        "sample_id": sample_id,
        "group": group,
        "total_fragments": total_fragments,
        "mapq30_fragments": mapq30_fragments,
        "usable_autosomal_fragments":
            usable_autosomal_fragments,
        "regional_fragments_100_220":
            regional_fragments,
        "mapq_threshold": MAPQ_MIN,
        "bin_size_bp": BIN_SIZE,
        "elapsed_minutes": elapsed_minutes,
    }])

    # Save only after successful completion.
    length_table.to_csv(
        global_output,
        index=False
    )

    regional_table.to_csv(
        regional_output,
        index=False
    )

    summary_table.to_csv(
        qc_output,
        index=False
    )

    print(
        f"Completed {sample_id}: "
        f"{total_fragments:,} fragments "
        f"in {elapsed_minutes:.2f} min"
    )


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    cancer_files = sorted(
        CANCER_DIR.glob("*.hg38.frag.tsv.bgz")
    )

    healthy_files = sorted(
        HEALTHY_DIR.glob("*.hg38.frag.tsv.bgz")
    )

    print("\nFiles detected:")
    print("Breast cancer:", len(cancer_files))
    print("Healthy:", len(healthy_files))

    # FULL COHORT PROCESSING

print("\nStarting full cohort processing...")

for i, file_path in enumerate(cancer_files, start=1):
    print(f"\nCancer sample {i}/{len(cancer_files)}")
    process_sample(file_path, "Breast cancer")

for i, file_path in enumerate(healthy_files, start=1):
    print(f"\nHealthy sample {i}/{len(healthy_files)}")
    process_sample(file_path, "Healthy")

print("\nFULL COHORT PROCESSING COMPLETE")