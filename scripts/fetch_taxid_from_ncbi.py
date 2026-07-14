#!/usr/bin/env python3
"""
Fetch correct TaxIDs from NCBI for all sequences in the database.
Uses NCBI E-utilities API for batch queries.
"""

import os
import sys
import time
import json
import urllib.request
import urllib.parse
from collections import defaultdict

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BATCH_SIZE = 200  # NCBI recommends max 200 IDs per request
DELAY = 0.4  # NCBI recommends 3 requests per second max

def fetch_taxid_batch(accessions, db="nucleotide"):
    """Fetch TaxIDs for a batch of accessions from NCBI."""
    ids = ",".join(accessions)
    
    # Step 1: esearch to get UIDs
    search_url = f"{NCBI_EUTILS_BASE}/esearch.fcgi?db={db}&term={ids}[Accn]&retmode=json&retmax={len(accessions)}"
    
    try:
        response = urllib.request.urlopen(search_url, timeout=30)
        data = json.loads(response.read().decode())
        uids = data.get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        print(f"Error in esearch: {e}")
        return {}
    
    if not uids:
        return {}
    
    # Step 2: efetch to get TaxIDs
    fetch_url = f"{NCBI_EUTILS_BASE}/efetch.fcgi?db={db}&id={','.join(uids)}&rettype=docsum&retmode=json"
    
    try:
        response = urllib.request.urlopen(fetch_url, timeout=60)
        data = json.loads(response.read().decode())
    except Exception as e:
        print(f"Error in efetch: {e}")
        return {}
    
    result = {}
    for item in data.get("result", {}).get("uids", []):
        acc = data["result"][item].get("accessionversion", "")
        taxid = data["result"][item].get("taxid", "")
        if acc and taxid:
            # Handle split sequences (e.g., GCA_xxx_1 -> GCA_xxx)
            base_acc = acc.split("_")[0] + "_" + acc.split("_")[1] if "_" in acc else acc
            result[acc] = taxid
            result[base_acc] = taxid
    
    return result

def fetch_assembly_taxid_batch(accessions):
    """Fetch TaxIDs for GenBank Assembly accessions (GCA_xxx) using batch queries."""
    result = {}
    batch_size = 50  # Smaller batch to avoid URL length limit
    
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]
        print(f"    Assembly batch {i//batch_size + 1}/{(len(accessions)//batch_size)+1}: {len(batch)} accessions")
        
        # Build search terms for this batch
        terms = [f"{a}[AssemblyAccn]" for a in batch]
        search_term = " OR ".join(terms)
        
        # Step 1: esearch to get UIDs
        search_url = f"{NCBI_EUTILS_BASE}/esearch.fcgi?db=assembly&term={urllib.parse.quote(search_term)}&retmode=json&retmax=100"
        
        try:
            response = urllib.request.urlopen(search_url, timeout=60)
            data = json.loads(response.read().decode())
            uids = data.get("esearchresult", {}).get("idlist", [])
        except Exception as e:
            print(f"      Error in assembly esearch: {e}")
            time.sleep(DELAY)
            continue
        
        if not uids:
            time.sleep(DELAY)
            continue
        
        # Step 2: esummary for this batch
        fetch_url = f"{NCBI_EUTILS_BASE}/esummary.fcgi?db=assembly&id={','.join(uids)}&retmode=json"
        
        try:
            response = urllib.request.urlopen(fetch_url, timeout=60)
            data = json.loads(response.read().decode())
            
            for uid in data.get("result", {}).get("uids", []):
                acc = data["result"][uid].get("assemblyacc", "")
                taxid = data["result"][uid].get("taxid", "")
                if acc and taxid:
                    result[acc] = taxid
                    result[get_base_accession(acc)] = taxid
        except Exception as e:
            print(f"      Error in assembly esummary: {e}")
        
        time.sleep(DELAY)
    
    return result

def get_base_accession(acc):
    """Get base accession for split sequences (e.g., GCA_xxx_1 -> GCA_xxx)."""
    parts = acc.split("_")
    if len(parts) >= 2:
        return parts[0] + "_" + parts[1]
    return acc

def main(accession_file, output_file):
    """Main function to process all accessions."""
    
    # Read accessions
    with open(accession_file, 'r') as f:
        accessions = [line.strip() for line in f]
    
    print(f"Total accessions: {len(accessions)}")
    
    # Separate by type
    gca_accessions = [a for a in accessions if a.startswith("GCA_") or a.startswith("GCF_")]
    other_accessions = [a for a in accessions if not (a.startswith("GCA_") or a.startswith("GCF_"))]
    
    print(f"GCA/GCF (Assembly) accessions: {len(gca_accessions)}")
    print(f"Other (Nucleotide) accessions: {len(other_accessions)}")
    
    taxid_map = {}
    
    # Process nucleotide accessions in batches
    print("\nFetching TaxIDs from Nucleotide database...")
    for i in range(0, len(other_accessions), BATCH_SIZE):
        batch = other_accessions[i:i+BATCH_SIZE]
        print(f"  Batch {i//BATCH_SIZE + 1}/{(len(other_accessions)//BATCH_SIZE)+1}: {len(batch)} accessions")
        
        result = fetch_taxid_batch(batch, db="nucleotide")
        taxid_map.update(result)
        
        time.sleep(DELAY)
    
    # Process assembly accessions (slower, one by one)
    print("\nFetching TaxIDs from Assembly database...")
    unique_gca = set()
    for a in gca_accessions:
        base = get_base_accession(a)
        unique_gca.add(base)
    
    print(f"  Unique GCA/GCF accessions: {len(unique_gca)}")
    
    gca_result = fetch_assembly_taxid_batch(list(unique_gca))
    taxid_map.update(gca_result)
    
    # Write output
    print(f"\nWriting {len(taxid_map)} mappings to {output_file}")
    with open(output_file, 'w') as f:
        for acc, taxid in taxid_map.items():
            f.write(f"{acc}\t{taxid}\n")
    
    # Report missing
    missing = [a for a in accessions if a not in taxid_map and get_base_accession(a) not in taxid_map]
    print(f"\nMissing TaxIDs: {len(missing)}")
    if missing[:10]:
        print(f"  First 10 missing: {missing[:10]}")
    
    return taxid_map

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python fetch_taxid_from_ncbi.py <accession_file> <output_file>")
        sys.exit(1)
    
    accession_file = sys.argv[1]
    output_file = sys.argv[2]
    
    main(accession_file, output_file)