"""VCF parsing and serialisation helpers.

These are the de-duplicated VCF utilities. In the original ``utils.py`` several
of these functions (notably :func:`process_vcf_dataframe`) were defined twice;
the *second*, working definition is the one preserved here. Clonality parsing
is factored into a standalone :func:`extract_clonality` so it can be tested in
isolation, and ``re`` is imported once at module scope rather than inside the
hot loop.
"""

from __future__ import annotations

import re
import warnings

import pandas as pd

__all__ = [
    "info_field_to_dict",
    "dict_to_vcf_info",
    "vcf_to_dataframe",
    "dataframe_to_vcf",
    "extract_clonality",
    "process_vcf_dataframe",
]


def _read_record(record) -> dict:
    """Flatten a :mod:`vcfpy` record into a plain dict."""
    return {
        "Chromosome": record.CHROM,
        "Position": record.POS,
        "ID": record.ID,
        "Reference Allele": record.REF,
        "Alternate Alleles": record.ALT,
        "Quality Score": record.QUAL,
        "Info": record.INFO,
    }


def info_field_to_dict(s: str) -> dict:
    """Parse a ``key:value;key:value`` INFO string into a typed dict."""
    d: dict = {}
    for item in s.split(";"):
        if not item:
            continue
        key, val = item.split(":", 1)
        if "," in val:
            d[key] = val.split(",")
        elif val == "TRUE":
            d[key] = True
        elif val == "FALSE":
            d[key] = False
        else:
            try:
                d[key] = int(val)
            except ValueError:
                try:
                    d[key] = float(val)
                except ValueError:
                    d[key] = val
    return d


def dict_to_vcf_info(d: dict) -> str:
    """Serialise an INFO dict back into the ``key:value;...`` string form."""
    parts = []
    for k, v in d.items():
        v_str = ",".join(str(i) for i in v) if isinstance(v, list) else str(v)
        parts.append(f"{k}:{v_str}")
    return ";".join(parts)


def vcf_to_dataframe(vcf_file: str) -> pd.DataFrame:
    """Read a VCF file into a DataFrame of raw records."""
    import vcfpy  # lazy import: keeps ``import hrdtimer`` lightweight

    # VCF contig headers often omit `length`; suppress the resulting vcfpy noise.
    # FieldInfoNotFound subclasses VCFPyWarning(Warning), not UserWarning.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=vcfpy.exceptions.FieldInfoNotFound)
        with open(vcf_file, "r") as f:
            reader = vcfpy.Reader.from_stream(f)
            data = [_read_record(record) for record in reader]
    return pd.DataFrame(data)


def dataframe_to_vcf(df: pd.DataFrame, output_vcf: str) -> None:
    """Write a processed mutation DataFrame back to a minimal VCF file."""
    cols = df.columns[: df.columns.get_loc("INFO") + 1]
    df = df[cols].copy()
    df["INFO"] = df["INFO"].apply(dict_to_vcf_info)
    with open(output_vcf, "w") as vcf:
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tINFO\n")
    df.to_csv(output_vcf, sep="\t", mode="a", index=False, header=False)


def extract_clonality(value) -> str | None:
    """Normalise a MutationTimeR ``CLS`` value to a clonality label.

    ``clonal [early]`` -> ``early``; ``subclonal`` -> ``subclonal``; anything
    else -> ``None`` (filtered downstream).
    """
    v = str(value)
    if v.startswith("clonal"):
        match = re.search(r"\[(.+?)\]", v)
        # Normalise to lowercase so "clonal [NA]" → "na" matches cls.lower() in
        # the preprocessor (MutationTimeR writes "NA" uppercase, others lowercase).
        return match.group(1).lower() if match else None
    if v == "subclonal":
        return "subclonal"
    return None


def process_vcf_dataframe(df: pd.DataFrame, time_analysis: bool = False) -> pd.DataFrame:
    """Extract timing-relevant fields from a raw VCF DataFrame.

    Parameters
    ----------
    df
        Output of :func:`vcf_to_dataframe`.
    time_analysis
        If ``True``, restrict to high-confidence timing regions
        (``MajCN == 2`` and ``inWGDregion``).
    """
    df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "INFO"]
    df["ID"] = df["ID"].apply(lambda x: x[0] if x else ".")
    df["ALT"] = df["ALT"].apply(lambda x: x[0].value if x else ".")

    fields = ["MajCN", "MinCN", "MutCN", "context", "CLS", "inWGDregion", "powr"]
    for field in fields:
        df[field] = df["INFO"].apply(lambda info: info.get(field, "."))

    df["CLS"] = df["CLS"].apply(extract_clonality)
    df = df[(df["context"] != ".") & df["CLS"].notna()]

    if time_analysis:
        df = df[(df["MajCN"] == 2) & (df["inWGDregion"] == "TRUE")]

    return df.rename(columns={"powr": "pow"})
