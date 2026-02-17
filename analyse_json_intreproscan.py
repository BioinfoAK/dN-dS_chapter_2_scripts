import json
import pandas as pd

def parse_interpro_with_cluster_ids(json_file, output_xlsx):
    with open(json_file, 'r') as f:
        data = json.load(f)

    extracted_rows = []

    for sequence in data.get('results', []):
        # --- NEW LOGIC TO FIND THE CLUSTER NAME ---
        seq_id = "unknown"
        xrefs = sequence.get('crossReferences', [])
        for ref in xrefs:
            ref_name = ref.get('name', '')
            if ref_name.startswith('cluster'):
                seq_id = ref_name
                break
        # ------------------------------------------
        
        orfs = sequence.get('openReadingFrames', [])
        if not orfs:
            continue
            
        sense_orfs = [o for o in orfs if o.get('strand') == 'SENSE']
        if not sense_orfs:
            continue
            
        longest_orf_obj = max(sense_orfs, key=lambda x: len(x.get('protein', {}).get('sequence', '')))
        protein_data = longest_orf_obj.get('protein', {})
        
        families = []
        domains = []
        
        for match in protein_data.get('matches', []):
            signature = match.get('signature', {})
            entry = signature.get('entry')
            
            if entry:
                entry_type = entry.get('type')
                row = {
                    'Nucleotide_Cluster_ID': seq_id,
                    'ORF_Length': len(protein_data.get('sequence', '')),
                    'Accession': entry.get('accession'),
                    'Description': entry.get('description'),
                    'Type': entry_type,
                    'Method': signature.get('signatureLibraryRelease', {}).get('library')
                }
                
                if entry_type == 'FAMILY':
                    families.append(row)
                elif entry_type == 'DOMAIN':
                    domains.append(row)

        final_list = families if families else domains
        extracted_rows.extend(final_list)

    df = pd.DataFrame(extracted_rows)
    
    if df.empty:
        print("Still nothing? Check if 'SENSE' is exactly as written in your JSON.")
    else:
        df.to_excel(output_xlsx, index=False)
        print(f"Slayed. Saved {len(df)} rows with Cluster IDs to {output_xlsx}")
        return df.head()

# Execute
df = parse_interpro_with_cluster_ids('iprscan5-R20260206-170212-0732-51481633-p1m.json', 'InterProScan_Summary_legprog.xlsx')
if df is not None:
    display(df)
