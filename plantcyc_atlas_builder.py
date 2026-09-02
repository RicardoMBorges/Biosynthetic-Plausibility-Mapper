#!/usr/bin/env python3
"""
plantcyc_atlas_builder.py

Convert a PlantCyc Pathway Tools flat-file archive into normalized tables and
curated reaction edges for the Biosynthetic Plausibility Mapper.

This version also creates a conservative chemical-identity layer so common
acid/base and salt forms can be reconciled across PlantCyc and LOTUS without
altering the original PlantCyc identifiers or structures.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import tarfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

try:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    RDKIT_AVAILABLE = True
except Exception:
    Chem = None
    rdMolStandardize = None
    RDKIT_AVAILABLE = False


FIELD_RE = re.compile(r"^([A-Z0-9?^_-]+) - (.*)$")
TAG_RE = re.compile(r"<[^>]+>")
FRAME_RE = re.compile(r"\|FRAME:\s*([^\s|]+)(?:\s+[^|]*)?\|")
CITS_RE = re.compile(r"\|CITS?:\s*([^|]+)\|")
TOKEN_RE = re.compile(r'"(?:\\.|[^"])*"|\|[^|]*\||[^\s()]+')


def clean_text(value: str) -> str:
    if value is None:
        return ""
    text = str(value)
    text = FRAME_RE.sub(lambda m: m.group(1), text)
    text = CITS_RE.sub("", text)
    text = TAG_RE.sub("", text)
    text = html.unescape(text)
    return re.sub(r"\s+", " ", text).strip()


def strip_token(token: str) -> str:
    token = token.strip()
    if len(token) >= 2 and token[0] == token[-1] == '"':
        token = token[1:-1]
    if len(token) >= 2 and token[0] == "|" and token[-1] == "|":
        token = token[1:-1]
    return token


def tokenize_atoms(text: str) -> list[str]:
    if not text:
        return []
    return [strip_token(t) for t in TOKEN_RE.findall(text)
            if strip_token(t) not in {"", "NIL"}]


def iter_dat_records(lines: Iterable[str]) -> Iterable[dict[str, list[str]]]:
    record: dict[str, list[str]] = defaultdict(list)
    last_field: str | None = None
    for raw_line in lines:
        line = raw_line.rstrip("\r\n")
        if not line or line.startswith("#"):
            continue
        if line == "//":
            if record:
                yield dict(record)
            record = defaultdict(list)
            last_field = None
            continue
        if line.startswith("/") and last_field and record.get(last_field):
            record[last_field][-1] = (record[last_field][-1] + " " + line[1:].strip()).strip()
            continue
        match = FIELD_RE.match(line)
        if match:
            field, value = match.groups()
            record[field].append(value)
            last_field = field
        elif last_field and record.get(last_field):
            record[last_field][-1] = (record[last_field][-1] + " " + line.strip()).strip()
    if record:
        yield dict(record)


def first(record: dict[str, list[str]], field: str, default: str = "") -> str:
    vals = record.get(field, [])
    return vals[0] if vals else default


def values(record: dict[str, list[str]], field: str) -> list[str]:
    return list(record.get(field, []))


REQUIRED_FILES = ["compounds.dat", "reactions.dat", "pathways.dat", "enzrxns.dat"]
OPTIONAL_FILES = ["proteins.dat", "genes.dat", "species.dat"]


def archive_member_map(tf: tarfile.TarFile) -> dict[str, str]:
    mapping = {}
    for member in tf.getmembers():
        basename = Path(member.name).name
        if basename in REQUIRED_FILES + OPTIONAL_FILES and "/data/" in member.name:
            mapping[basename] = member.name
    return mapping


def read_dat_from_tar(tf: tarfile.TarFile, member_name: str):
    extracted = tf.extractfile(member_name)
    if extracted is None:
        raise FileNotFoundError(member_name)
    text = extracted.read().decode("utf-8", errors="replace")
    return list(iter_dat_records(text.splitlines()))


def parse_formula(parts: list[str]) -> str:
    atoms = []
    for part in parts:
        m = re.match(r"^\(\s*([A-Za-z][A-Za-z0-9]*)\s+([0-9]+)\s*\)$", part)
        if m:
            element, count = m.groups()
            atoms.append((element, count))
    return "".join(element + (count if count != "1" else "") for element, count in atoms)


def parse_dblinks(entries: list[str]) -> dict[str, list[str]]:
    output = defaultdict(list)
    for entry in entries:
        toks = tokenize_atoms(entry)
        if len(toks) >= 2:
            db, accession = toks[0], toks[1]
            if accession and accession not in output[db]:
                output[db].append(accession)
    return dict(output)


def serialize_list(items: Iterable[Any]) -> str:
    return "|".join(str(x) for x in items if str(x).strip())


def serialize_dblinks(dblinks: dict[str, list[str]]) -> str:
    chunks = []
    for db in sorted(dblinks):
        for accession in dblinks[db]:
            chunks.append(f"{db}:{accession}")
    return "|".join(chunks)


def parse_reaction_layout(layout: str) -> dict[str, Any] | None:
    if not layout:
        return None
    first_match = re.match(r'^\(\s*("[^"]+"|\|[^|]+\||[^\s()]+)', layout)
    if not first_match:
        return None
    reaction_id = strip_token(first_match.group(1))
    left_match = re.search(r"\(:LEFT-PRIMARIES\s+([^)]*)\)", layout)
    right_match = re.search(r"\(:RIGHT-PRIMARIES\s+([^)]*)\)", layout)
    direction_match = re.search(r"\(:DIRECTION\s+:([A-Z0-9-]+)\)", layout)
    return {
        "reaction_id": reaction_id,
        "left_primaries": tokenize_atoms(left_match.group(1)) if left_match else [],
        "right_primaries": tokenize_atoms(right_match.group(1)) if right_match else [],
        "direction": direction_match.group(1) if direction_match else "",
        "raw_layout": layout,
    }


def normalize_reaction_direction(value: str) -> str:
    v = str(value or "").upper()
    if v in {"PHYSIOL-LEFT-TO-RIGHT", "LEFT-TO-RIGHT", "L2R"}:
        return "L2R"
    if v in {"PHYSIOL-RIGHT-TO-LEFT", "RIGHT-TO-LEFT", "R2L"}:
        return "R2L"
    if "REVERS" in v:
        return "REVERSIBLE"
    return ""


def direction_to_pairs(left, right, direction):
    d = normalize_reaction_direction(direction)
    if d == "R2L":
        return right, left, True
    if d == "L2R":
        return left, right, True
    return left, right, False


# ---------------------------------------------------------------------
# Chemical identity normalization
# ---------------------------------------------------------------------

def valid_standard_inchikey(value: Any) -> str:
    raw = str(value or "").strip().upper()
    if raw.startswith("INCHIKEY="):
        raw = raw.split("=", 1)[1]
    parts = raw.split("-")
    if len(parts) == 3 and len(parts[0]) == 14 and len(parts[1]) == 10 and len(parts[2]) == 1:
        return raw
    return ""


def inchikey_connectivity_key(value: Any) -> str:
    key = valid_standard_inchikey(value)
    return key.split("-", 1)[0] if key else ""


def inchikey_protonation_key(value: Any) -> str:
    key = valid_standard_inchikey(value)
    if not key:
        return ""
    parts = key.split("-")
    return f"{parts[0]}-{parts[1]}"


def normalize_identity_structure(smiles: Any) -> dict[str, str]:
    result = {
        "identity_parent_smiles": "",
        "identity_parent_inchi": "",
        "identity_parent_inchikey": "",
        "identity_parent_connectivity_key": "",
        "identity_parent_protonation_key": "",
        "identity_parent_smiles_no_stereo": "",
        "identity_normalization_status": "not_available",
    }

    raw = str(smiles or "").strip()
    if not raw:
        result["identity_normalization_status"] = "missing_smiles"
        return result
    if not RDKIT_AVAILABLE:
        result["identity_normalization_status"] = "rdkit_unavailable"
        return result
    if any(token in raw for token in ("[R]", "[R1]", "[R2]", "[R3]", "[R4]")):
        result["identity_normalization_status"] = "generic_structure"
        return result

    try:
        mol = Chem.MolFromSmiles(raw, sanitize=True)
    except Exception:
        mol = None
    if mol is None:
        result["identity_normalization_status"] = "smiles_parse_failed"
        return result

    try:
        parent = rdMolStandardize.FragmentParent(mol)
        parent = rdMolStandardize.Uncharger().uncharge(parent)

        # Rebuild from SMILES so RDKit caches/ring information are initialized.
        tmp = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        parent = Chem.MolFromSmiles(tmp, sanitize=True)
        if parent is None:
            raise ValueError("Could not rebuild normalized parent")
        parent.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(parent)

        parent_smiles = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        parent_inchi = Chem.MolToInchi(parent)
        parent_inchikey = Chem.InchiToInchiKey(parent_inchi) if parent_inchi else ""

        no_stereo = Chem.Mol(parent)
        Chem.RemoveStereochemistry(no_stereo)
        parent_smiles_no_stereo = Chem.MolToSmiles(
            no_stereo, canonical=True, isomericSmiles=False
        )

        result.update({
            "identity_parent_smiles": parent_smiles,
            "identity_parent_inchi": parent_inchi,
            "identity_parent_inchikey": parent_inchikey,
            "identity_parent_connectivity_key": inchikey_connectivity_key(parent_inchikey),
            "identity_parent_protonation_key": inchikey_protonation_key(parent_inchikey),
            "identity_parent_smiles_no_stereo": parent_smiles_no_stereo,
            "identity_normalization_status": "normalized",
        })
    except Exception:
        result["identity_normalization_status"] = "parent_normalization_failed"

    return result


def build_compounds(records):
    rows, lookup = [], {}
    for rec in records:
        compound_id = first(rec, "UNIQUE-ID")
        if not compound_id:
            continue
        dblinks = parse_dblinks(values(rec, "DBLINKS"))
        name_raw = first(rec, "COMMON-NAME", compound_id)
        inchi_key = first(rec, "INCHI-KEY")
        if inchi_key.startswith("InChIKey="):
            inchi_key = inchi_key.split("=", 1)[1]

        original_smiles = first(rec, "SMILES")
        identity = normalize_identity_structure(original_smiles)

        row = {
            "compound_id": compound_id,
            "name": clean_text(name_raw) or compound_id,
            "common_name_raw": name_raw,
            "synonyms": serialize_list(clean_text(x) for x in values(rec, "SYNONYMS")),
            "types": serialize_list(values(rec, "TYPES")),
            "formula": parse_formula(values(rec, "CHEMICAL-FORMULA")),
            "molecular_weight": first(rec, "MOLECULAR-WEIGHT"),
            "smiles": original_smiles,
            "inchi": first(rec, "INCHI"),
            "non_standard_inchi": first(rec, "NON-STANDARD-INCHI"),
            "inchikey": inchi_key,
            "source_inchikey_connectivity_key": inchikey_connectivity_key(inchi_key),
            "source_inchikey_protonation_key": inchikey_protonation_key(inchi_key),
            **identity,
            "chebi_id": "|".join(dblinks.get("CHEBI", [])),
            "kegg_compound_id": "|".join(dblinks.get("LIGAND-CPD", [])),
            "metanetx_id": "|".join(dblinks.get("METANETX", [])),
            "dblinks": serialize_dblinks(dblinks),
            "source_database": "PlantCyc",
        }
        rows.append(row)
        lookup[compound_id] = row
    return rows, lookup


def build_proteins(records):
    rows, lookup = [], {}
    for rec in records:
        protein_id = first(rec, "UNIQUE-ID")
        if not protein_id:
            continue
        name = clean_text(first(rec, "COMMON-NAME"))
        if not name:
            syn = values(rec, "SYNONYMS")
            name = clean_text(syn[0]) if syn else protein_id
        row = {
            "protein_id": protein_id,
            "name": name,
            "gene_ids": serialize_list(values(rec, "GENE")),
            "catalyzes": serialize_list(values(rec, "CATALYZES")),
            "species_ids": serialize_list(values(rec, "SPECIES")),
            "synonyms": serialize_list(clean_text(x) for x in values(rec, "SYNONYMS")),
            "source_database": "PlantCyc",
        }
        rows.append(row)
        lookup[protein_id] = row
    return rows, lookup


def build_enzymatic_reactions(records, protein_lookup=None):
    rows = []
    by_reaction = defaultdict(list)
    protein_lookup = protein_lookup or {}
    for rec in records:
        enzrxn_id = first(rec, "UNIQUE-ID")
        reaction_ids = values(rec, "REACTION")
        enzyme_ids = values(rec, "ENZYME")
        enzyme_names = []
        for enzyme_id in enzyme_ids:
            protein = protein_lookup.get(enzyme_id, {})
            enzyme_names.append(
                protein.get("name") or clean_text(first(rec, "COMMON-NAME")) or enzyme_id
            )
        row = {
            "enzrxn_id": enzrxn_id,
            "name": clean_text(first(rec, "COMMON-NAME")),
            "reaction_ids": serialize_list(reaction_ids),
            "enzyme_ids": serialize_list(enzyme_ids),
            "enzyme_names": serialize_list(enzyme_names),
            "reaction_direction": first(rec, "REACTION-DIRECTION"),
            "basis_for_assignment": first(rec, "BASIS-FOR-ASSIGNMENT"),
            "physiologically_relevant": first(rec, "PHYSIOLOGICALLY-RELEVANT?"),
            "citations": serialize_list(values(rec, "CITATIONS")),
            "source_database": "PlantCyc",
        }
        rows.append(row)
        for reaction_id in reaction_ids:
            by_reaction[reaction_id].append(row)
    return rows, by_reaction


def build_reactions(records, enzrxns_by_reaction):
    rows, lookup = [], {}
    for rec in records:
        reaction_id = first(rec, "UNIQUE-ID")
        if not reaction_id:
            continue
        enz_rows = enzrxns_by_reaction.get(reaction_id, [])
        enzyme_names, enzrxn_ids = [], []
        for er in enz_rows:
            enzrxn_ids.append(er["enzrxn_id"])
            enzyme_names.extend(x for x in str(er["enzyme_names"]).split("|") if x)
        ec_numbers = [x.removeprefix("EC-") for x in values(rec, "EC-NUMBER")]
        row = {
            "reaction_id": reaction_id,
            "name": clean_text(first(rec, "COMMON-NAME")) or reaction_id,
            "types": serialize_list(values(rec, "TYPES")),
            "left_ids": serialize_list(values(rec, "LEFT")),
            "right_ids": serialize_list(values(rec, "RIGHT")),
            "reaction_direction": first(rec, "REACTION-DIRECTION"),
            "pathway_ids": serialize_list(values(rec, "IN-PATHWAY")),
            "enzymatic_reaction_ids": serialize_list(enzrxn_ids),
            "enzyme_names": serialize_list(dict.fromkeys(enzyme_names)),
            "ec_numbers": serialize_list(ec_numbers),
            "citations": serialize_list(values(rec, "CITATIONS")),
            "physiologically_relevant": first(rec, "PHYSIOLOGICALLY-RELEVANT?"),
            "orphan": first(rec, "ORPHAN?"),
            "reaction_balance_status": first(rec, "REACTION-BALANCE-STATUS"),
            "source_database": "PlantCyc",
        }
        rows.append(row)
        lookup[reaction_id] = row
    return rows, lookup


def build_pathways(records):
    pathway_rows, lookup, pathway_reaction_rows = [], {}, []
    for rec in records:
        pathway_id = first(rec, "UNIQUE-ID")
        if not pathway_id:
            continue
        reaction_list = values(rec, "REACTION-LIST")
        layouts = [p for x in values(rec, "REACTION-LAYOUT")
                   if (p := parse_reaction_layout(x))]
        row = {
            "pathway_id": pathway_id,
            "name": clean_text(first(rec, "COMMON-NAME")) or pathway_id,
            "types": serialize_list(values(rec, "TYPES")),
            "reaction_ids": serialize_list(reaction_list),
            "species_ids": serialize_list(values(rec, "SPECIES")),
            "taxonomic_range": serialize_list(values(rec, "TAXONOMIC-RANGE")),
            "super_pathways": serialize_list(values(rec, "SUPER-PATHWAYS")),
            "sub_pathways": serialize_list(values(rec, "SUB-PATHWAYS")),
            "citations": serialize_list(values(rec, "CITATIONS")),
            "dblinks": serialize_dblinks(parse_dblinks(values(rec, "DBLINKS"))),
            "reaction_layout_count": len(layouts),
            "source_database": "PlantCyc",
        }
        pathway_rows.append(row)
        lookup[pathway_id] = row

        layout_by_rxn = {x["reaction_id"]: x for x in layouts}
        all_rxns = list(dict.fromkeys(reaction_list + list(layout_by_rxn)))
        for position, reaction_id in enumerate(all_rxns, start=1):
            layout = layout_by_rxn.get(reaction_id)
            pathway_reaction_rows.append({
                "pathway_id": pathway_id,
                "pathway_name": row["name"],
                "reaction_id": reaction_id,
                "position_hint": position,
                "has_reaction_layout": bool(layout),
                "layout_direction": layout["direction"] if layout else "",
                "left_primaries": serialize_list(layout["left_primaries"]) if layout else "",
                "right_primaries": serialize_list(layout["right_primaries"]) if layout else "",
                "raw_layout": layout["raw_layout"] if layout else "",
                "source_database": "PlantCyc",
            })
    return pathway_rows, lookup, pathway_reaction_rows


def compound_label(compound_id, compounds):
    c = compounds.get(compound_id)
    return (c.get("name") if c else "") or compound_id


def reaction_enzyme_label(reaction_id, reactions):
    return str(reactions.get(reaction_id, {}).get("enzyme_names", ""))


def reaction_name(reaction_id, reactions):
    r = reactions.get(reaction_id, {})
    return str(r.get("name", "")) or reaction_id


def edge_row(pathway_id, pathway_name, reaction_id, source_id, target_id,
             directed, basis, evidence, confidence, reactions, compounds, direction):
    sc = compounds.get(source_id, {})
    tc = compounds.get(target_id, {})
    return {
        "source": compound_label(source_id, compounds),
        "target": compound_label(target_id, compounds),
        "reaction": reaction_name(reaction_id, reactions),
        "enzyme": reaction_enzyme_label(reaction_id, reactions),
        "pathway": pathway_name,
        "evidence": evidence,
        "confidence": confidence,
        "directed": directed,
        "source_id": source_id,
        "target_id": target_id,
        "source_inchikey": sc.get("inchikey", ""),
        "target_inchikey": tc.get("inchikey", ""),
        "source_smiles": sc.get("smiles", ""),
        "target_smiles": tc.get("smiles", ""),
        "source_formula": sc.get("formula", ""),
        "target_formula": tc.get("formula", ""),
        "source_identity_parent_smiles": sc.get("identity_parent_smiles", ""),
        "target_identity_parent_smiles": tc.get("identity_parent_smiles", ""),
        "source_identity_parent_inchikey": sc.get("identity_parent_inchikey", ""),
        "target_identity_parent_inchikey": tc.get("identity_parent_inchikey", ""),
        "source_identity_connectivity_key": sc.get("identity_parent_connectivity_key", ""),
        "target_identity_connectivity_key": tc.get("identity_parent_connectivity_key", ""),
        "source_identity_protonation_key": sc.get("identity_parent_protonation_key", ""),
        "target_identity_protonation_key": tc.get("identity_parent_protonation_key", ""),
        "reaction_id": reaction_id,
        "ec_numbers": reactions.get(reaction_id, {}).get("ec_numbers", ""),
        "pathway_id": pathway_id,
        "direction": normalize_reaction_direction(direction),
        "edge_basis": basis,
        "citations": reactions.get(reaction_id, {}).get("citations", ""),
        "source_database": "PlantCyc",
    }


def build_edges(pathway_records, pathways, reactions, compounds,
                fallback_to_reaction_sides=True, pathway_filter=None):
    edges, unresolved = [], []
    seen = set()

    for rec in pathway_records:
        pathway_id = first(rec, "UNIQUE-ID")
        if not pathway_id:
            continue
        if pathway_filter and pathway_id != pathway_filter:
            continue

        pathway_name = pathways.get(pathway_id, {}).get("name", pathway_id)
        layout_entries = [p for x in values(rec, "REACTION-LAYOUT")
                          if (p := parse_reaction_layout(x))]
        layout_rxn_ids = {p["reaction_id"] for p in layout_entries}

        for layout in layout_entries:
            reaction_id = layout["reaction_id"]
            sources, targets, directed = direction_to_pairs(
                layout["left_primaries"], layout["right_primaries"], layout["direction"]
            )
            if not sources or not targets:
                unresolved.append({
                    "pathway_id": pathway_id,
                    "reaction_id": reaction_id,
                    "reason": "REACTION-LAYOUT missing primary source or target",
                    "raw": layout["raw_layout"],
                })
                continue
            for source_id in sources:
                for target_id in targets:
                    key = (pathway_id, reaction_id, source_id, target_id, "reaction-layout")
                    if key in seen:
                        continue
                    seen.add(key)
                    edges.append(edge_row(
                        pathway_id, pathway_name, reaction_id, source_id, target_id,
                        directed, "reaction-layout", "PlantCyc REACTION-LAYOUT",
                        1.0, reactions, compounds, layout["direction"]
                    ))

        if not fallback_to_reaction_sides:
            continue

        for reaction_id in values(rec, "REACTION-LIST"):
            if reaction_id in layout_rxn_ids:
                continue
            reaction = reactions.get(reaction_id)
            if not reaction:
                unresolved.append({
                    "pathway_id": pathway_id,
                    "reaction_id": reaction_id,
                    "reason": "reaction not found in reactions.dat",
                    "raw": "",
                })
                continue
            left = [x for x in str(reaction.get("left_ids", "")).split("|") if x]
            right = [x for x in str(reaction.get("right_ids", "")).split("|") if x]
            sources, targets, directed = direction_to_pairs(
                left, right, str(reaction.get("reaction_direction", ""))
            )
            for source_id in sources:
                for target_id in targets:
                    key = (pathway_id, reaction_id, source_id, target_id, "reaction-sides-fallback")
                    if key in seen:
                        continue
                    seen.add(key)
                    edges.append(edge_row(
                        pathway_id, pathway_name, reaction_id, source_id, target_id,
                        directed, "reaction-sides-fallback", "PlantCyc LEFT/RIGHT fallback",
                        0.75, reactions, compounds, reaction.get("reaction_direction", "")
                    ))

    return edges, unresolved


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    columns, seen = [], set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                columns.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_database(archive: Path, output_dir: Path,
                   fallback_to_reaction_sides=True, pathway_filter=None):
    output_dir.mkdir(parents=True, exist_ok=True)

    with tarfile.open(archive, "r:gz") as tf:
        members = archive_member_map(tf)
        missing = [x for x in REQUIRED_FILES if x not in members]
        if missing:
            raise FileNotFoundError(
                "Required PlantCyc files missing from archive: " + ", ".join(missing)
            )
        compound_records = read_dat_from_tar(tf, members["compounds.dat"])
        reaction_records = read_dat_from_tar(tf, members["reactions.dat"])
        pathway_records = read_dat_from_tar(tf, members["pathways.dat"])
        enzrxn_records = read_dat_from_tar(tf, members["enzrxns.dat"])
        protein_records = (
            read_dat_from_tar(tf, members["proteins.dat"])
            if "proteins.dat" in members else []
        )

    protein_rows, protein_lookup = build_proteins(protein_records)
    compound_rows, compound_lookup = build_compounds(compound_records)
    enzrxn_rows, enzrxns_by_reaction = build_enzymatic_reactions(
        enzrxn_records, protein_lookup
    )
    reaction_rows, reaction_lookup = build_reactions(
        reaction_records, enzrxns_by_reaction
    )
    pathway_rows, pathway_lookup, pathway_reaction_rows = build_pathways(pathway_records)
    edge_rows, unresolved_rows = build_edges(
        pathway_records, pathway_lookup, reaction_lookup, compound_lookup,
        fallback_to_reaction_sides=fallback_to_reaction_sides,
        pathway_filter=pathway_filter
    )

    write_csv(output_dir / "atlas_compounds.csv", compound_rows)
    write_csv(output_dir / "atlas_reactions.csv", reaction_rows)
    write_csv(output_dir / "atlas_pathways.csv", pathway_rows)
    write_csv(output_dir / "atlas_pathway_reactions.csv", pathway_reaction_rows)
    write_csv(output_dir / "atlas_enzymatic_reactions.csv", enzrxn_rows)
    if protein_rows:
        write_csv(output_dir / "atlas_proteins.csv", protein_rows)
    write_csv(output_dir / "atlas_curated_edges.csv", edge_rows)
    write_csv(output_dir / "atlas_unresolved_edges.csv", unresolved_rows)

    report = {
        "archive": str(archive),
        "output_dir": str(output_dir),
        "rdkit_available": RDKIT_AVAILABLE,
        "pathway_filter": pathway_filter,
        "fallback_to_reaction_sides": fallback_to_reaction_sides,
        "counts": {
            "compounds": len(compound_rows),
            "identity_normalized_compounds": sum(
                r.get("identity_normalization_status") == "normalized"
                for r in compound_rows
            ),
            "reactions": len(reaction_rows),
            "pathways": len(pathway_rows),
            "pathway_reactions": len(pathway_reaction_rows),
            "enzymatic_reactions": len(enzrxn_rows),
            "proteins": len(protein_rows),
            "curated_edges": len(edge_rows),
            "unresolved_edges": len(unresolved_rows),
            "edges_from_reaction_layout": sum(
                r.get("edge_basis") == "reaction-layout" for r in edge_rows
            ),
            "edges_from_fallback": sum(
                r.get("edge_basis") == "reaction-sides-fallback" for r in edge_rows
            ),
        }
    }
    (output_dir / "build_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return report


def cli():
    parser = argparse.ArgumentParser(
        description="Build Biosynthetic Plausibility Mapper tables from PlantCyc tar.gz."
    )
    parser.add_argument("archive", type=Path)
    parser.add_argument("-o", "--output-dir", type=Path,
                        default=Path("plantcyc_atlas_db"))
    parser.add_argument("--no-fallback", action="store_true")
    parser.add_argument("--pathway", default=None)
    return parser.parse_args()


def main():
    args = cli()
    if not args.archive.exists():
        raise SystemExit(f"Archive not found: {args.archive}")

    report = build_database(
        args.archive,
        args.output_dir,
        fallback_to_reaction_sides=not args.no_fallback,
        pathway_filter=args.pathway,
    )

    print("\nPlantCyc -> Biosynthetic Plausibility Mapper build complete")
    print("=" * 62)
    print(f"RDKit identity normalization: {'enabled' if RDKIT_AVAILABLE else 'disabled'}")
    for key, value in report["counts"].items():
        print(f"{key:32s}: {value:,}")
    print(f"\nOutput directory: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
