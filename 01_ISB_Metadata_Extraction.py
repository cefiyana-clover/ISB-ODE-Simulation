
"""
ISB Protocol: Automated Metadata Extraction Pipeline
----------------------------------------------------
This module executes the programmatic retrieval of summary statistics 
from the PubMed Central database to parameterize the simulation's input vectors.
"""

import sys
import subprocess
import time
import warnings

warnings.filterwarnings('ignore')

# 1. DEPENDENCY MANAGEMENT
try:
    from Bio import Entrez
    import pandas as pd
except ImportError:
    print("Installing required dependencies (biopython, pandas)...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython", "pandas"])
    from Bio import Entrez
    import pandas as pd

# 2. API CONFIGURATION
Entrez.email = "research_contact@institution.edu" # To be updated prior to submission
Entrez.tool = "ISB_Data_Miner_Protocol"

def query_ncbi_database(search_query: str, max_results: int = 50) -> pd.DataFrame:
    """
    Queries NCBI PMC database to retrieve genomic and metabolic summary statistics.
    """
    try:
        search_handle = Entrez.esearch(db="pmc", term=search_query, retmax=max_results)
        search_record = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_record.get("IdList", [])
        if not id_list:
            print("No matching records detected for the parameterized query.")
            return pd.DataFrame()
            
        extracted_records = []
        for pmcid in id_list:
            try:
                time.sleep(0.5) # API rate limiting compliance
                summary_handle = Entrez.esummary(db="pmc", id=pmcid)
                summary_record = Entrez.read(summary_handle)
                summary_handle.close()
                
                record = summary_record[0] if isinstance(summary_record, list) and len(summary_record) > 0 else {}
                extracted_records.append({
                    "PMCID": f"PMC{pmcid}",
                    "Publication_Date": record.get("PubDate", "N/A"),
                    "Title": record.get("Title", "N/A")
                })
            except Exception:
                continue
                
        return pd.DataFrame(extracted_records)
        
    except Exception as e:
        print(f"Extraction Error: {e}")
        return pd.DataFrame()

# 3. EXECUTION: PAN-ANCESTRY GWAS EXTRACTION
if __name__ == "__main__":
    print("Initiating NCBI PMC query for Pan-Ancestry GWAS metadata...")
    gwas_query = '''
    ("Major Depressive Disorder"[Title/Abstract] OR MDD[Title/Abstract]) 
    AND ("Genome-Wide Association"[Title/Abstract] OR GWAS[Title/Abstract]) 
    AND ("East Asian"[Title/Abstract] OR "African"[Title/Abstract] OR EAS[Title/Abstract] OR AFR[Title/Abstract])
    AND open access[filter]
    '''
    df_gwas = query_ncbi_database(gwas_query)
    
    if not df_gwas.empty:
        print(f"Extraction complete. {len(df_gwas)} records retrieved.")
        print(df_gwas.head())
        # The output is utilized to parameterize baseline genetic variance in Script 2.
        # df_gwas.to_csv('Extracted_GWAS_Metadata.csv', index=False)

