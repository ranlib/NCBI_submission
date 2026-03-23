#!/usr/bin/env python3
# ==============================================================================
# Script Name: Wastewater Metadata Test (one row per sample-target)
# Purpose:
#   1) Read RunName from sequencing SampleSheet.csv
#   2) Make a WBE-only copy of the sample sheet named: <RunName>-SampleSheet.csv 
#   3) Query the CHHS API once per requested target
#   4) Write one output row per WBE sample per requested target specified
#   5) Save all generated CSVs in a local Metadata-Output f older
#
# ==============================================================================

import csv
import json
import re
from datetime import datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen


# Input SampleSheet in the current working directory.
SAMPLESHEET_PATH = Path("SampleSheet.csv")

# Folder where all generated CSV files will be written.
OUTPUT_DIR = Path("Metadata-Output")

# CHHS CKAN API settings.
API_URL = "https://data.chhs.ca.gov/api/3/action/datastore_search"
RESOURCE_ID = "2742b824-3736-4292-90a9-7fad98e94c06"
PAGE_LIMIT = 1000

# Prefix added to the output sample_name column.
SAMPLE_NAME_PREFIX = "DWRL_"

# Base filters applied to every API request.
DB_BASE_FILTERS = {
    "lab_id": "DWRL",
}

# ------------------------------------------------------------------------------
# TARGET CONFIG
# Add or remove targets here.
# Example:
# REQUESTED_TARGETS = [
#     {"pcr_target": "sars-cov-2"},
#     {"pcr_target": "rsv"}, #add targets as listed in database
# ]
# ------------------------------------------------------------------------------
REQUESTED_TARGETS = [
    {"pcr_target": "sars-cov-2"},
    #{"pcr_target": "rsv"}, 
]


# ------------------------------------------------------------------------------
# OUTPUT COLUMNS
# Comment out any line you do not want in the output CSVs.
#
# source_type options:
#   sample_name          -> raw WBE sample name from SampleSheet.csv
#   db                   -> field from the matched database row
#   fixed                -> hard-coded text
#   formula              -> custom function in FORMULAS
#   target_request       -> field from REQUESTED_TARGETS entry
#   db_or_target_request -> use db field if present, otherwise requested target field
# ------------------------------------------------------------------------------
OUTPUT_COLUMNS = [
    ("sample_name", "formula", "prefixed_sample_name"), #Required
    ("ncbi-spuid", "db", "sample_id"), #Required
    ("ncbi-spuid-sra", "sample_name", ""), #Required
    ("ncbi-bioproject", "fixed", "PRJNA747181"), #Required
    ("isolation_source", "db", "sample_matrix"), #Required
    ("organism", "fixed", "wastewater metagenomics"), #Required
    ("collection_date", "db", "sample_collect_date"), #Required
    ("collection_time", "db", "sample_collect_time"), 
    ("country", "db", "county_treatmentplant"), #Required
    ("state", "db", "reporting_jurisdiction"), #Required
    ("collection_site_id", "db", "site_id"),
    ("ww_sample_site", "db", "site_id"),
    ("ww_flow", "db", "flow_rate"),
    ("instantaneous_flow", "db", "flow_rate"),
    ("ww_population", "db", "population_served"), #Required
    ("ww_surv_jurisdiction", "db", "reporting_jurisdiction"),
    ("ww_sample_matrix", "db", "sample_matrix"), #Required
    ("ww_sample_type", "db", "sample_type"), #Required
    ("ww_sample_duration", "formula", "ww_sample_duration"), #Required
    ("ww_temperature", "db", "collection_water_temp"),
    ("ww_ph", "db", "ph"),
    ("ww_industrial_effluent_percent", "db", "industrial_input"),
    ("ww_surv_system_sample_id", "db", "sample_id"),
    ("ww_pre_treatment", "db", "pretreatment"),
    ("specimen_processing", "db", "major_lab_method"),
    ("specimen_processing_details", "db", "major_lab_method_desc"),
    ("ww_processing_protocol", "db", "major_lab_method_desc"),
    ("concentration_method", "db", "concentration_method"),
    ("extraction_method", "db", "extraction_method"),
    ("extraction_control", "db", "rec_eff_target_name"),
    ("ww_endog_control_1", "db", "hum_frac_target_mic"),
    ("ww_endog_control_1_conc", "db", "hum_frac_mic_conc"),
    ("ww_endog_control_1_units", "db", "hum_frac_mic_unit"),
    ("ww_endog_control_2", "db", "hum_frac_target_chem"),
    ("ww_endog_control_2_conc", "db", "hum_frac_chem_conc"),
    ("ww_endog_control_2_units", "db", "hum_frac_chem_unit"),
    ("ww_surv_target_1", "db_or_target_request", ("pcr_target", "pcr_target")), #Required
    ("ww_surv_target_1_known_present", "formula", "target_known_present"), #Required
    ("ww_surv_target_1_protocol", "db", "validated_methods"),
    ("ww_surv_target_1_conc", "db", "pcr_target_avg_conc"),
    ("ww_surv_target_1_conc_unit", "db", "pcr_target_units"),
    ("ww_surv_target_1_gene", "db", "pcr_gene_target"),
    ("sequenced_by", "fixed", "CDPH-CLS-Genomics Center"),
]


# ------------------------------------------------------------------------------
# FORMULAS
# These helper functions fill derived output values.
# ------------------------------------------------------------------------------
def formula_prefixed_sample_name(sample_name, db_row, target_request):
    """Add the configured prefix to the raw WBE sample name."""
    if not sample_name:
        return ""
    return f"{SAMPLE_NAME_PREFIX}{sample_name}"


def formula_ww_sample_duration(sample_name, db_row, target_request):
    """Return 24 for composite samples and 0 for grab samples."""
    if not db_row:
        return ""

    sample_type = str(db_row.get("sample_type", "")).strip().lower()
    if sample_type == "composite" or "composite" in sample_type:
        return "24"
    if sample_type == "grab" or "grab" in sample_type:
        return "0"
    return ""


def formula_target_known_present(sample_name, db_row, target_request):
    """Mark the target as present only when a database match exists."""
    return "Yes" if db_row else ""


FORMULAS = {
    "prefixed_sample_name": formula_prefixed_sample_name,
    "ww_sample_duration": formula_ww_sample_duration,
    "target_known_present": formula_target_known_present,
}


# ------------------------------------------------------------------------------
# SampleSheet helpers
# These functions read the Illumina SampleSheet and keep only WBE rows.
# ------------------------------------------------------------------------------
def is_blank_row(row):
    """Return True when a CSV row is empty or all blank cells."""
    return (not row) or all(cell.strip() == "" for cell in row)


def read_csv_rows(path):
    """Read a CSV-like file into a list of raw rows."""
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.reader(handle))


def get_run_name(rows):
    """Extract RunName from the SampleSheet header block."""
    for row in rows:
        if row and row[0].strip() == "RunName" and len(row) > 1:
            run_name = row[1].strip()
            if run_name:
                return run_name
    raise ValueError("Could not find RunName in SampleSheet.csv")


def find_section(rows, section_name):
    """Return the row index for a section header like [BCLConvert_Data]."""
    for i, row in enumerate(rows):
        if row and row[0].strip() == section_name:
            return i
    return None


def get_section_table(rows, section_name):
    """Return header and row bounds for a section table."""
    section_idx = find_section(rows, section_name)
    if section_idx is None:
        return None

    header_idx = section_idx + 1
    while header_idx < len(rows) and is_blank_row(rows[header_idx]):
        header_idx += 1

    if header_idx >= len(rows):
        raise ValueError(f"Could not find header row under {section_name}")

    header = [cell.strip() for cell in rows[header_idx]]
    if "Sample_ID" not in header:
        raise ValueError(f"{section_name} is missing the Sample_ID column")

    data_start_idx = header_idx + 1
    data_end_idx = data_start_idx
    while data_end_idx < len(rows):
        row = rows[data_end_idx]
        if row and row[0].strip().startswith("["):
            break
        data_end_idx += 1

    return header, data_start_idx, data_end_idx


def filter_section_to_wbe(rows, section_name, keep_names=None):
    """Keep only WBE rows, or only rows whose Sample_ID is in keep_names."""
    section = get_section_table(rows, section_name)
    if section is None:
        return rows, []

    header, data_start_idx, data_end_idx = section
    sample_id_idx = header.index("Sample_ID")

    kept_rows = []
    kept_names_in_order = []

    for i in range(data_start_idx, data_end_idx):
        row = rows[i]
        if is_blank_row(row):
            continue

        padded = [cell.strip() for cell in row] + [""] * max(0, len(header) - len(row))
        padded = padded[: len(header)]
        sample_name = padded[sample_id_idx].strip()

        if keep_names is None:
            keep = "WBE" in sample_name.upper()
        else:
            keep = sample_name in keep_names

        if keep:
            kept_rows.append(padded)
            kept_names_in_order.append(sample_name)

    new_rows = []
    new_rows.extend(rows[:data_start_idx])
    new_rows.extend(kept_rows)
    new_rows.extend(rows[data_end_idx:])
    return new_rows, kept_names_in_order


def unique_in_order(items):
    """Remove duplicates while preserving the original order."""
    seen = set()
    out = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


# ------------------------------------------------------------------------------
# Database helpers
# These functions fetch and normalize CHHS database records.
# ------------------------------------------------------------------------------
def normalize_text(value):
    """Trim and lowercase a value for stable comparisons."""
    return str(value or "").strip().lower()


def normalize_db_sample_name(raw_sample_id):
    """Keep only the part of sample_id after the last underscore."""
    raw_sample_id = str(raw_sample_id or "").strip()
    return re.sub(r"^.*_", "", raw_sample_id).strip()


def normalize_date(value):
    """Convert known date formats to YYYY-MM-DD when possible."""
    value = str(value or "").strip()
    if not value:
        return ""

    for fmt in (
        "%Y-%m-%d",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%dT%H:%M:%S.%f",
        "%m/%d/%Y",
        "%m/%d/%y",
    ):
        try:
            return datetime.strptime(value.replace("Z", ""), fmt).date().isoformat()
        except ValueError:
            pass

    try:
        return datetime.fromisoformat(value.replace("Z", "")).date().isoformat()
    except Exception:
        return value


def normalize_time(value):
    """Convert known time formats to HH:MM:SS when possible."""
    value = str(value or "").strip()
    if not value:
        return ""

    for fmt in (
        "%H:%M:%S",
        "%H:%M",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%dT%H:%M:%S.%f",
    ):
        try:
            return datetime.strptime(value.replace("Z", ""), fmt).time().isoformat()
        except ValueError:
            pass

    try:
        return datetime.fromisoformat(value.replace("Z", "")).time().isoformat()
    except Exception:
        return value


def clean_db_record(record):
    """Normalize one raw API record for downstream matching and export."""
    cleaned = {}
    for key, value in record.items():
        cleaned[key] = "" if value is None else str(value).strip()

    cleaned["sample_name"] = normalize_db_sample_name(cleaned.get("sample_id", ""))
    cleaned["sample_collect_date"] = normalize_date(cleaned.get("sample_collect_date", ""))
    cleaned["sample_collect_time"] = normalize_time(cleaned.get("sample_collect_time", ""))
    cleaned["sample_type_norm"] = normalize_text(cleaned.get("sample_type", ""))
    cleaned["pcr_target_norm"] = normalize_text(cleaned.get("pcr_target", ""))

    if cleaned.get("sample_type") in (
        "24-hr time-weighted composite",
        "24-hr flow-weighted composite",
    ):
        cleaned["sample_type"] = "composite"
        cleaned["sample_type_norm"] = "composite"

    return cleaned


def normalize_requested_targets(targets):
    """Add lowercase comparison keys to each requested target."""
    normalized = []
    for target in targets:
        normalized.append({
            "pcr_target": str(target.get("pcr_target", "")).strip(),
            "pcr_target_norm": normalize_text(target.get("pcr_target", "")),
        })
    return normalized


NORMALIZED_REQUESTED_TARGETS = normalize_requested_targets(REQUESTED_TARGETS)


def fetch_target_records(target_request):
    """Fetch all API records for one requested pcr_target."""
    offset = 0
    all_records = []

    filters = dict(DB_BASE_FILTERS)
    filters["pcr_target"] = target_request["pcr_target"]

    while True:
        params = {
            "resource_id": RESOURCE_ID,
            "filters": json.dumps(filters),
            "limit": PAGE_LIMIT,
            "offset": offset,
        }
        url = f"{API_URL}?{urlencode(params)}"

        with urlopen(url) as response:
            payload = json.loads(response.read().decode("utf-8"))

        if not payload.get("success", False):
            raise RuntimeError(f"API request failed: {payload}")

        result = payload.get("result", {})
        records = result.get("records", [])
        total = int(result.get("total", 0))

        if not records:
            break

        all_records.extend(records)
        offset += len(records)

        if offset >= total:
            break

    return all_records


def build_db_index_by_sample_and_target():
    """
    Build a nested lookup like:
        {sample_name: {pcr_target_norm: [cleaned_row, ...]}}
    """
    index = {}

    for target_request in NORMALIZED_REQUESTED_TARGETS:
        records = fetch_target_records(target_request)
        for record in records:
            cleaned = clean_db_record(record)
            sample_name = cleaned.get("sample_name", "")
            pcr_target_norm = cleaned.get("pcr_target_norm", "")

            if not sample_name or not pcr_target_norm:
                continue

            sample_bucket = index.setdefault(sample_name, {})
            sample_bucket.setdefault(pcr_target_norm, []).append(cleaned)

    return index


def choose_target_row(candidate_rows):
    """Pick the first matching row for a sample + target."""
    if not candidate_rows:
        return None
    return candidate_rows[0]


# ------------------------------------------------------------------------------
# Output helpers
# These functions map database fields into the requested CSV columns.
# ------------------------------------------------------------------------------
def value_from_spec(sample_name, db_row, target_request, spec):
    """Resolve one output value from the configured column spec."""
    _output_header, source_type, source_value = spec

    if source_type == "sample_name":
        return sample_name

    if source_type == "db":
        if not db_row:
            return ""
        return str(db_row.get(source_value, "")).strip()

    if source_type == "fixed":
        return str(source_value)

    if source_type == "formula":
        func = FORMULAS.get(source_value)
        if not func:
            return ""
        return str(func(sample_name, db_row, target_request) or "")

    if source_type == "target_request":
        return str(target_request.get(source_value, "")).strip()

    if source_type == "db_or_target_request":
        db_field, request_field = source_value
        if db_row and str(db_row.get(db_field, "")).strip():
            return str(db_row.get(db_field, "")).strip()
        return str(target_request.get(request_field, "")).strip()

    raise ValueError(f"Unknown source_type: {source_type}")


def build_output_row(sample_name, db_row, target_request, column_specs):
    """Build one output row using OUTPUT_COLUMNS."""
    row = {}
    for spec in column_specs:
        header = spec[0]
        row[header] = value_from_spec(sample_name, db_row, target_request, spec)
    return row


def write_output_csv(path, rows, column_specs):
    """Write a list of dictionaries to CSV using configured headers."""
    fieldnames = [spec[0] for spec in column_specs]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
def main():
    """Run the full WBE filtering, API matching, and CSV export workflow."""
    if not SAMPLESHEET_PATH.exists():
        raise FileNotFoundError(f"Could not find {SAMPLESHEET_PATH.resolve()}")

    # Create the output folder if it does not already exist.
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Read the input SampleSheet and find the run name.
    original_rows = read_csv_rows(SAMPLESHEET_PATH)
    run_name = get_run_name(original_rows)

    # Keep only WBE rows in the BCLConvert section.
    filtered_rows, wbe_names = filter_section_to_wbe(
        original_rows,
        "[BCLConvert_Data]",
        keep_names=None,
    )
    wbe_names = unique_in_order(wbe_names)

    # Keep the same WBE rows in Cloud_Data, if that section exists.
    filtered_rows, _ = filter_section_to_wbe(
        filtered_rows,
        "[Cloud_Data]",
        keep_names=set(wbe_names),
    )

    # Write the WBE-only SampleSheet copy.
    output_samplesheet = OUTPUT_DIR / f"{run_name}-SampleSheet.csv"
    with output_samplesheet.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(filtered_rows)

    # Fetch database rows and index them by sample name and target.
    db_index = build_db_index_by_sample_and_target()

    reportable_rows = []
    not_reported_rows = []

    # Create one output row per sample per requested target.
    for sample_name in wbe_names:
        sample_targets = db_index.get(sample_name, {})

        for target_request in NORMALIZED_REQUESTED_TARGETS:
            candidate_rows = sample_targets.get(target_request["pcr_target_norm"], [])
            matched_row = choose_target_row(candidate_rows)

            if matched_row:
                reportable_rows.append(
                    build_output_row(sample_name, matched_row, target_request, OUTPUT_COLUMNS)
                )
            else:
                not_reported_rows.append(
                    build_output_row(sample_name, None, target_request, OUTPUT_COLUMNS)
                )

    # Write both report CSVs to Metadata-Output.
    reportable_path = OUTPUT_DIR / f"reportable_{run_name}.csv"
    not_reported_path = OUTPUT_DIR / f"not-reported_{run_name}.csv"

    write_output_csv(reportable_path, reportable_rows, OUTPUT_COLUMNS)
    write_output_csv(not_reported_path, not_reported_rows, OUTPUT_COLUMNS)

    print(f"RunName:            {run_name}")
    print(f"WBE samples found:  {len(wbe_names)}")
    print(f"Output folder:      {OUTPUT_DIR.resolve()}")
    print(f"Copied SampleSheet: {output_samplesheet.name}")
    print(f"Reportable CSV:     {reportable_path.name} ({len(reportable_rows)} rows)")
    print(f"Not reported CSV:   {not_reported_path.name} ({len(not_reported_rows)} rows)")
    print("Requested targets:  " + ", ".join(t["pcr_target"] for t in REQUESTED_TARGETS))


if __name__ == "__main__":
    try:
        main()
    except HTTPError as error:
        print(f"HTTP error: {error.code} {error.reason}")
    except URLError as error:
        print(f"URL error: {error.reason}")
    except Exception as error:
        print(f"Unexpected error: {error}")
