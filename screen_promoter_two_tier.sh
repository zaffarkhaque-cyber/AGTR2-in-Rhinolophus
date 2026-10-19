python
#!/usr/bin/env python
import os
import sys
import hashlib
import zipfile
import shutil
import tempfile
import subprocess

# ----------------- CONSTANTS & EMBEDDED QUERIES -----------------
QUERY_MYOTIS_SEQ = (
    "GGTTTAAGCACTTTAAATACATTAGCTTGTTTAATTTTCCTGACTCTATGATGTAGGTAA"
    "TTAACATCATCCCCATTTATCAGATGAGAAAGCTGAGACACAAAAGACGGAAATAACTTG"
    "GCCAAAGTCATAGAAGTAGTAAGTGGTATAGCTGAGCTTGTGAAACCAGGCAGTTGGATC"
    "CAAGGTCATTTACATACCAGTATGCTACATTGCCTCTCTATTAAATTACCGAGCATAAAG"
    "CATCTAGCACAGTGTTTGGCACATAGTTTTTCAACTAGTATTATCACTACCATCACCTTT"
    "TCTTTCTTCCTTCCTTTTGTTACTCTTGCTTCTCTTGGGCTGATCTTGTAAGGTTGATCC"
    "CTGCAGCGCCCACAAGAAGAGATAACAATTCTATCTAGGCTGCAGTTGCTGCCTTTCTCC"
    "TTCTTCACCTTGTGCACTTTGTATTTAAACCGCAAATATGACAGAGACCTATGAAAGGGA"
    "ACACCAAATGTTAATACCTCTTTTCCAAAACATGAAAATTATATCCCTCCTATCACCCCC"
    "GCTCCAGCAACAAGCCAACAGAAATTGTTTAAAGCAAGGGAGAAATTTTAGAAATTCTAA"
    "ACTTCAAGTCACATTGCTACTGTTCTTTAATAAACAAGGGTAAAGATTTTTAAAAGTAAT"
    "GAAGGTTACAAAGAACATTGGGGGACAATGACACTGCTACTAGAATAAATTCCCTAGTCT"
    "CCACTCCTTAATGCGTTTTAGAATGCCAGTATCAGGGAAGAAGAGAAAATTCTACAAATA"
    "ATTATATTCTTCTTTTCAAGTAAGTCACTGACTTCAAAATTCTCTAACACATTTGTGGAA"
    "ACTTCATTTTTTTTGTTTGAGATTTATTTGAATGAGCTGTTATGATTGGAGACAGTGAGA"
    "ATTTCAGATTAATGTTTTGCAGCCAAAAAAAAAAGCTCTGGAAAGCTGGCAAGTGTTCAT"
    "AAGTCAGCCCTAGAATTATGTAGGTTGAAGGCTCCCAGTGG"
)

QUERY_APUS_SEQ = (
    "TCTCTTCTGGTTCACCTCACTGTAAGCAGAGCTGGTTATGGAGCCTCTGATGCCTGTTGC"
    "CTCTGGTGTGCTACCAGGTCCAAGAGCAAGGGCAGGCACCAGGATTCCCCCAAGCCACCA"
    "CAACTGGATACCCACCCACACTGCAAGGCTGCAGTGTTCAAAGGTGTTTAAAGAATGGGA"
    "TGAGCTTGAGTCTGGCCCTGGCTGGTGATGAGCTGGGGCAAAGGGACATGTCACTGACAG"
    "GGACTTCTGTGGTTGCTTTCCACGTGACATGTGCAGCAATGCCTGGAGAGCACAGCACAC"
    "ATGTCCACAGTAGATGGGCCACTGCTGCAGAGGACATCACGCTCATGATATCAGGGCAAA"
    "ATATTGAAGGAGATGCTTATTTTGGCCTTGTGTTTATCCCCCAGCAGCCCTCTTGCAGCT"
    "TCTTCCATGGAGCCGAGACTGGGGCCAGGTGGACCACCAGGAAGAAGCAGATGTAAGAAC"
    "TCAGTTCTGTGTGAGTGGGAAGTAGCTGTGAGCCTGCAGCCTCGTCACCTGCTTGCTTTG"
    "GTCCTGTGCTACCCCAGTGCAGACTCAGAGCCTCTGCTCTGCCACGTGGATTTGTCCAGC"
    "TGCCCCCCTCAGTTCCAGCAGTGTCCCCCCACAGGGACTTGAGCAGTGGCTTTCCCCTCC"
    "TGCTGATGGGGAGTCAGGATCAGCCCCCTCACACAGGCAATACTGGTGTGAATGCTGGGA"
    "TGTTCCTCCAGTCTCTGTCTGTGTGTGTTCACGTTGTGGTCAGTCCTGCCAGGACTAACA"
    "GCCTATTCTCCTCTCAAATGTGTTAGCAAGTGCTGGGGCTCCAGGGAGTGCATGGGAATA"
    "CCCCTGCCCTCAGGGATCCATGGATCCTTTCACAAAGGATCACCTGGGGGACAGAAGAGG"
    "GTTGCATGTGGGTGTATCCTCCTGGCCCAAGAGCTGCAGGAGAAAACCCACAGAGATATT"
    "TGTGTGGATTTCCTTTTTCTTTCTCTCTTTTTTTTTTTCT"
)

# ----------------- INPUT DATA CONFIGURATION (WSL Mapped) -----------------
INPUT_FILES = [
    "/mnt/c/Users/shazi/Downloads/Rhinolophus_ferrumequinum_ncbi_dataset.zip",
    "/mnt/c/Users/shazi/Downloads/Rhinolophus_hipposideros_ncbi_dataset (1).zip",
    "/mnt/c/Users/shazi/Downloads/Rhinolophus_affinis_ncbi_dataset.zip",
    "/mnt/c/Users/shazi/Downloads/Rhinolophus_sinicus_ncbi_dataset.zip",
    "/mnt/c/Users/shazi/OneDrive/Desktop/BATCH GENOMES/Group 5/Hipposideros_larvatus_ncbi_dataset.zip",
    "/mnt/c/Users/shazi/Downloads/Hipposideros_armiger_ncbi_dataset.zip",
    "/mnt/c/Users/shazi/Downloads/Hipposideros_larvatus_ncbi_dataset.zip",
    "/mnt/c/Users/shazi/OneDrive/Desktop/BATCH GENOMES/Group7/Myotis_myotis_ncbi_dataset.zip",
    "/mnt/c/Users/shazi/Downloads/Calypte_anna_ncbi_dataset (2).zip",
    "/mnt/c/Users/shazi/Downloads/Apus_apus_ncbi_dataset.zip",
    "/mnt/c/Users/shazi/Downloads/Selasphorus_platycercus_ncbi_dataset.zip"
]

OUTPUT_TSV = "promoter_screen_output.tsv"

def get_md5(file_path):
    hasher = hashlib.md5()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b''):
            hasher.update(chunk)
    return hasher.hexdigest()

def clean_species_name(filename):
    name = filename.replace("_ncbi_dataset", "").replace(".zip", "")
    if " (" in name:
        name = name.split(" (")[0]
    return name.strip().replace(" ", "_")

def classify_hit(alignment_len, pident, qlen):
    if not alignment_len or alignment_len == 0:
        return "NO_HIT"
    elif (alignment_len * 100 >= 60 * qlen) and (pident >= 80.0):
        return "RETAINED_CANDIDATE"
    elif (alignment_len * 100 >= 30 * qlen):
        return "PARTIAL_FRAGMENT"
    else:
        return "SHORT_REMNANT"

def extract_largest_fasta_from_zip(zip_path, extract_to):
    largest_file = None
    largest_size = 0
    with zipfile.ZipFile(zip_path, 'r') as archive:
        for info in archive.infolist():
            if info.filename.endswith(('.fna', '.fa', '.fasta')) and info.file_size > largest_size:
                largest_size = info.file_size
                largest_file = info.filename
        if largest_file:
            base_name = os.path.basename(largest_file)
            dest_path = os.path.join(extract_to, base_name)
            with archive.open(largest_file) as source, open(dest_path, 'wb') as target:
                shutil.copyfileobj(source, target)
            return dest_path
    return None

def run_blast_pass(query_path, db_path, task="blastn"):
    cmd = [
        "blastn", "-query", query_path, "-db", db_path,
        "-task", task, "-outfmt", "6 sseqid length pident evalue sstart send",
        "-max_target_seqs", "1"
    ]
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        lines = res.stdout.strip().split('\n')
        if lines and lines[0]:
            parts = lines[0].split('\t')
            return {
                "sseqid": parts[0], "length": int(parts[1]), "pident": float(parts[2]),
                "evalue": parts[3], "region": f"{parts[4]}-{parts[5]}"
            }
    except subprocess.CalledProcessError:
        pass
    return None

def main():
    print("Starting Two-Tiered Promoter Screening Pipeline...")
    seen_hashes = {}
    
    with tempfile.TemporaryDirectory() as scratch_dir:
        q_myotis_path = os.path.join(scratch_dir, "query_myotis.fa")
        q_apus_path = os.path.join(scratch_dir, "query_apus.fa")
        
        with open(q_myotis_path, 'w') as f:
            f.write(f">MYOTIS_MYOTIS_AGTR2_PROMOTER\n{QUERY_MYOTIS_SEQ}\n")
        with open(q_apus_path, 'w') as f:
            f.write(f">APUS_APUS_AGTR2_PROMOTER\n{QUERY_APUS_SEQ}\n")
            
        len_myotis = len(QUERY_MYOTIS_SEQ)
        len_apus = len(QUERY_APUS_SEQ)
        
        queries = [
            {"name": "Myotis promoter", "path": q_myotis_path, "len": len_myotis},
            {"name": "Apus promoter", "path": q_apus_path, "len": len_apus}
        ]

        with open(OUTPUT_TSV, 'w') as out_tsv:
            out_tsv.write("Species\tSource file\tQuery\tBest length (bp)\tMax length (bp)\t% ident\tE-value\tSubject region\tClassification\n")

            for file_path in INPUT_FILES:
                if not os.path.exists(file_path):
                    print(f"⚠️ Warning: File path does not exist, skipping: {file_path}")
                    continue
                
                raw_filename = os.path.basename(file_path)
                file_hash = get_md5(file_path)
                
                if file_hash in seen_hashes:
                    print(f"🛑 SKIP (byte-identical duplicate of {seen_hashes[file_hash]}): {raw_filename}")
                    continue
                seen_hashes[file_hash] = raw_filename
                
                species_token = clean_species_name(raw_filename)
                print(f"🧬 Processing: {species_token} ({raw_filename})")
                
                genome_fa = extract_largest_fasta_from_zip(file_path, scratch_dir)
                if not genome_fa:
                    print(f"❌ No valid FASTA files found in {raw_filename}")
                    continue
                
                db_prefix = os.path.join(scratch_dir, f"db_{species_token}")
                make_db_cmd = ["makeblastdb", "-in", genome_fa, "-dbtype", "nucl", "-out", db_prefix]
                subprocess.run(make_db_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
                
                for q in queries:
                    hit = run_blast_pass(q["path"], db_prefix, task="blastn")
                    if not hit:
                        hit = run_blast_pass(q["path"], db_prefix, task="blastn-short")
                        
                    if hit:
                        classification = classify_hit(hit["length"], hit["pident"], q["len"])
                        out_tsv.write(f"{species_token}\t{raw_filename}\t{q['name']}\t{hit['length']}\t{hit['length']}\t{hit['pident']}\t{hit['evalue']}\t{hit['region']}\t{classification}\n")
                    else:
                        out_tsv.write(f"{species_token}\t{raw_filename}\t{q['name']}\t0\t0\t0.0\tN/A\tN/A\tNO_HIT\n")
                
                try: os.remove(genome_fa)
                except: pass

    print(f"\n🎉 Screening Finished! Results written to: '{os.path.abspath(OUTPUT_TSV)}'")

if __name__ == "__main__":
    main()




    To extract the sequence from Hipposideros larvatus  that matched:

    
    #!/usr/bin/env python
import os
import zipfile
import subprocess

zip_path = "/mnt/c/Users/shazi/Downloads/Hipposideros_larvatus_ncbi_dataset.zip"
query_seq = (
    "GGTTTAAGCACTTTAAATACATTAGCTTGTTTAATTTTCCTGACTCTATGATGTAGGTAA"
    "TTAACATCATCCCCATTTATCAGATGAGAAAGCTGAGACACAAAAGACGGAAATAACTTG"
    "GCCAAAGTCATAGAAGTAGTAAGTGGTATAGCTGAGCTTGTGAAACCAGGCAGTTGGATC"
    "CAAGGTCATTTACATACCAGTATGCTACATTGCCTCTCTATTAAATTACCGAGCATAAAG"
    "CATCTAGCACAGTGTTTGGCACATAGTTTTTCAACTAGTATTATCACTACCATCACCTTT"
    "TCTTTCTTCCTTCCTTTTGTTACTCTTGCTTCTCTTGGGCTGATCTTGTAAGGTTGATCC"
    "CTGCAGCGCCCACAAGAAGAGATAACAATTCTATCTAGGCTGCAGTTGCTGCCTTTCTCC"
    "TTCTTCACCTTGTGCACTTTGTATTTAAACCGCAAATATGACAGAGACCTATGAAAGGGA"
    "ACACCAAATGTTAATACCTCTTTTCCAAAACATGAAAATTATATCCCTCCTATCACCCCC"
    "GCTCCAGCAACAAGCCAACAGAAATTGTTTAAAGCAAGGGAGAAATTTTAGAAATTCTAA"
    "ACTTCAAGTCACATTGCTACTGTTCTTTAATAAACAAGGGTAAAGATTTTTAAAAGTAAT"
    "GAAGGTTACAAAGAACATTGGGGGACAATGACACTGCTACTAGAATAAATTCCCTAGTCT"
    "CCACTCCTTAATGCGTTTTAGAATGCCAGTATCAGGGAAGAAGAGAAAATTCTACAAATA"
    "ATTATATTCTTCTTTTCAAGTAAGTCACTGACTTCAAAATTCTCTAACACATTTGTGGAA"
    "ACTTCATTTTTTTTGTTTGAGATTTATTTGAATGAGCTGTTATGATTGGAGACAGTGAGA"
    "ATTTCAGATTAATGTTTTGCAGCCAAAAAAAAAAGCTCTGGAAAGCTGGCAAGTGTTCAT"
    "AAGTCAGCCCTAGAATTATGTAGGTTGAAGGCTCCCAGTGG"
)

# Write the query out temporarily
with open("temp_query.fa", "w") as f:
    f.write(f">MYOTIS_MYOTIS_AGTR2_PROMOTER\n{query_seq}\n")

print("Extracting genome database file...")
with zipfile.ZipFile(zip_path, 'r') as archive:
    largest_file = max([i for i in archive.infolist() if i.filename.endswith(('.fna', '.fa', '.fasta'))], key=lambda x: x.file_size)
    with archive.open(largest_file.filename) as src, open("temp_genome.fna", "wb") as dst:
        dst.write(src.read())

print("Building temporary local database...")
subprocess.run(["makeblastdb", "-in", "temp_genome.fna", "-dbtype", "nucl", "-out", "temp_db"], stdout=subprocess.DEVNULL)

# We ask BLAST to pull the matching subject sequence directly (sseq)
print("\n=== HIPPOSIDEROS MATCHING RETAINED PROMOTER SEQUENCE ===")
cmd = ["blastn", "-query", "temp_query.fa", "-db", "temp_db", "-outfmt", "6 sseqid sstart send sseq", "-max_target_seqs", "1"]
res = subprocess.run(cmd, text=True, capture_output=True)

if res.stdout.strip():
    parts = res.stdout.strip().split('\t')
    print(f"Chromosome/Scaffold ID: {parts[0]}")
    print(f"Genomic Coordinates: {parts[1]} to {parts[2]}")
    print("\n>Hipposideros_larvatus_retained_promoter")
    print(parts[3])
else:
    print("No matching alignment sequence found.")

# Clean up local workspace footprints
for f in ["temp_query.fa", "temp_genome.fna", "temp_db.nhr", "temp_db.nin", "temp_db.nsq"]:
    if os.path.exists(f): os.remove(f)
