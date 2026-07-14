#!/usr/bin/env python3
"""Test NCBI API for PP729064"""

import urllib.request
import json

acc = 'PP729064'
search_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term={acc}[Accn]&retmode=json'
response = urllib.request.urlopen(search_url, timeout=30)
data = json.loads(response.read().decode())
uids = data.get('esearchresult', {}).get('idlist', [])
print(f'PP729064 UID: {uids}')

if uids:
    fetch_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={uids[0]}&rettype=docsum&retmode=json'
    response = urllib.request.urlopen(fetch_url, timeout=30)
    data = json.loads(response.read().decode())
    for item in data.get('result', {}).get('uids', []):
        accver = data['result'][item].get('accessionversion', '')
        taxid = data['result'][item].get('taxid', '')
        print(f'Accession: {accver}, TaxID: {taxid}')