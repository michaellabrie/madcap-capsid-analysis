import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    from pathlib import Path

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    return Path, np, pd


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Sequence-Based Discovery and Prioritization of Novel AAV Capsids

    ## Research Question

    How does sequence similarity to known natural AAV isolates vary among
    MADCAP-discovered capsids, and can sequence-diversity metrics be used
    to identify candidates for further structural and experimental
    characterization?
    """)
    return


@app.cell
def _(Path, pd):
    data_path = Path(
        "data/raw/MADCAP_Table_S5_identity_matrix.csv"
    )

    identity_df = pd.read_csv(
        data_path,
        header=1
    )

    identity_df = identity_df.rename(
        columns={"Unnamed: 0": "capsid"}
    )

    print("Rows:", identity_df.shape[0])
    print("Columns:", identity_df.shape[1])

    identity_df.head()
    return (identity_df,)


@app.cell
def _(identity_df):
    print("First 20 capsids:")
    print(identity_df["capsid"].head(20).tolist())

    print("\nLast 20 capsids:")
    print(identity_df["capsid"].tail(20).tolist())
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Identifying Novel and Known AAV Capsids

    Capsid names beginning with "MC" represent capsids discovered through
    the MADCAP pipeline. The remaining capsids represent previously known
    natural AAV isolates.

    These groups are identified from their sequence names so that novel
    MADCAP capsids can be compared directly with known AAVs.
    """)
    return


@app.cell
def _(identity_df):
    # Clean any extra spaces from capsid names
    identity_df.columns = [
        str(column).strip()
        for column in identity_df.columns
    ]

    identity_df["capsid"] = (
        identity_df["capsid"]
        .astype(str)
        .str.strip()
    )

    # Separate MADCAP and previously known capsids
    madcap_capsids = [
        capsid
        for capsid in identity_df["capsid"]
        if capsid.startswith("MC")
    ]

    known_capsids = [
        capsid
        for capsid in identity_df["capsid"]
        if not capsid.startswith("MC")
    ]

    print("MADCAP capsids:", len(madcap_capsids))
    print("Known capsids:", len(known_capsids))
    print("Total capsids:", len(madcap_capsids) + len(known_capsids))
    return known_capsids, madcap_capsids


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Measuring Similarity to Known AAV Capsids

    For each MADCAP capsid, sequence identity is compared against all
    previously known capsids in the dataset.

    Three summary measurements are calculated:

    - nearest known capsid
    - highest sequence identity to a known capsid
    - average sequence identity to all known capsids

    Sequence divergence from the nearest known capsid is also calculated
    as 100 minus the maximum sequence identity.
    """)
    return


@app.cell
def _(identity_df, pd):
    # Keep only rows corresponding to MADCAP capsids
    madcap_df = identity_df[
        identity_df["capsid"].str.startswith("MC")
    ].copy()

    # Columns representing known AAV capsids
    known_columns = [
        column
        for column in identity_df.columns
        if column != "capsid"
        and not column.startswith("MC")
    ]

    # Make sure sequence-identity values are numeric
    known_identity_matrix = madcap_df[
        known_columns
    ].apply(
        pd.to_numeric,
        errors="coerce"
    )

    print("MADCAP capsids analyzed:", known_identity_matrix.shape[0])
    print("Known capsids compared:", known_identity_matrix.shape[1])
    return known_identity_matrix, madcap_df


@app.cell
def _(known_identity_matrix, madcap_df):
    madcap_analysis = madcap_df[["capsid"]].copy()

    madcap_analysis["nearest_known_capsid"] = (
        known_identity_matrix.idxmax(axis=1)
    )

    madcap_analysis["nearest_known_identity"] = (
        known_identity_matrix.max(axis=1)
    )

    madcap_analysis["mean_known_identity"] = (
        known_identity_matrix.mean(axis=1)
    )

    madcap_analysis["median_known_identity"] = (
        known_identity_matrix.median(axis=1)
    )

    madcap_analysis["sequence_divergence"] = (
        100 - madcap_analysis["nearest_known_identity"]
    )

    madcap_analysis.head()
    return (madcap_analysis,)


@app.cell
def _(madcap_analysis):
    most_divergent = (
        madcap_analysis
        .sort_values(
            "sequence_divergence",
            ascending=False
        )
    )

    most_divergent.head(15)
    return (most_divergent,)


@app.cell
def _(madcap_analysis, np):
    # Classify capsids using the divergence threshold used in the MADCAP study
    madcap_analysis["divergence_status"] = np.where(
        madcap_analysis["nearest_known_identity"] < 95,
        "Divergent (<95% identity)",
        "Similar to known isolate (≥95% identity)"
    )

    madcap_analysis["divergence_status"].value_counts()
    return


@app.cell
def _(madcap_analysis):
    # Summarize sequence-similarity measurements by divergence status
    divergence_summary = (
        madcap_analysis
        .groupby(
            "divergence_status",
            as_index=False
        )
        .agg(
            number_of_capsids=("capsid", "count"),
            mean_nearest_identity=("nearest_known_identity", "mean"),
            mean_overall_identity=("mean_known_identity", "mean"),
            mean_sequence_divergence=("sequence_divergence", "mean")
        )
    )

    divergence_summary
    return


@app.cell
def _(known_capsids, madcap_capsids):
    print("MADCAP capsids:", len(madcap_capsids))
    print("Known capsids:", len(known_capsids))
    return


@app.cell
def _(most_divergent):
    most_divergent.head(15)
    return


@app.cell
def _(madcap_analysis):
    # Check whether our calculated identity range agrees with the published study
    print(
        "Lowest nearest-known identity:",
        madcap_analysis["nearest_known_identity"].min()
    )

    print(
        "Highest nearest-known identity:",
        madcap_analysis["nearest_known_identity"].max()
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Validation Against the Published Study

    The calculated nearest-known sequence identities range from approximately
    58.4% to 99.9%. This agrees with the range reported by Lyashenko et al.,
    supporting that the sequence identity matrix has been interpreted correctly.
    """)
    return


if __name__ == "__main__":
    app.run()
