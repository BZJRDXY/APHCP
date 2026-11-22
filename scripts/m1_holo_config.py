"""Read the same query window and cutoff for each M1 stage."""

from datetime import date
import json
import hashlib


RELEASE_FIELD = "rcsb_accession_info.initial_release_date"


def read_config(run):
    return json.loads((run / "inputs" / "query.json").read_text(encoding="utf-8"))


def window(config):
    filters = date_filters(config["m1_1_request"])
    start = date.fromisoformat(filters["greater_or_equal"]["value"])
    end = date.fromisoformat(filters["less"]["value"])
    if start >= end:
        raise ValueError("The M1 window start must precede its exclusive end")
    cutoff_value = config.get("cutoff_date")
    return start, end, date.fromisoformat(cutoff_value) if cutoff_value else None


def date_filters(request):
    filters = {}
    for node in request["query"]["nodes"]:
        parameters = node.get("parameters", {})
        if parameters.get("attribute") == RELEASE_FIELD:
            operator = parameters["operator"]
            if operator in filters:
                raise ValueError(f"Duplicate M1 release date filter: {operator}")
            filters[operator] = parameters
    if set(filters) != {"greater_or_equal", "less"}:
        raise ValueError("M1 requires one inclusive start and one exclusive end filter")
    return filters



def read_json(path):
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temp.replace(path)


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def search_complete(search):
    ids = search.get("entry_ids", [])
    return (search.get("concordant") is True
            and search.get("total_count") == search.get("unique_ids") == len(ids) == len(set(ids)))


def metadata_complete(search, data):
    unavailable = data.get("confirmed_unavailable", {})
    return (search_complete(search) and not data.get("missing_ids", [None])
            and data.get("requested") == len(search["entry_ids"])
            and data.get("received", -1) + len(unavailable) == data["requested"]
            and set(unavailable) <= set(search["entry_ids"]))


def reusable_window(run):
    search = read_json(run / "search" / "summary.json") or {}
    receipt = read_json(run / "run_receipt.json") or {}
    if not search_complete(search):
        return False
    if not search["entry_ids"]:
        return True  # Legacy empty windows have only a concordant Search summary.
    data = read_json(run / "data" / "summary.json") or {}
    return (receipt.get("status") == "INVENTORY_COMPLETE" and metadata_complete(search, data)
            and (run / "inventory" / "candidate_entries.jsonl").is_file())
