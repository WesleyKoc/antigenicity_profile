import csv
import random
from datetime import datetime, timedelta

def parse_fasta(fasta_path):
    """Parses a FASTA file into a list of (header_name, sequence) tuples."""
    records = []
    current_name = None
    current_seq = []
    
    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_name:
                    records.append((current_name, "".join(current_seq)))
                current_name = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            records.append((current_name, "".join(current_seq)))
            
    return records

def extract_year_from_name(name):
    """Extracts a 4-digit year from the FASTA header name."""
    parts = name.split('/')
    for part in reversed(parts):
        if part.isdigit() and len(part) == 4:
            return int(part)
    return 2020  # Fallback year if none found

def generate_random_date(year):
    """Generates a random YYYY-MM-DD date for a given year."""
    start_date = datetime(year, 1, 1)
    end_date = datetime(year, 12, 31)
    delta = end_date - start_date
    random_days = random.randrange(delta.days + 1)
    random_date = start_date + timedelta(days=random_days)
    return random_date.strftime("%Y-%m-%d")

def generate_epi():
    """Generates a random EPI accession number."""
    return f"EPI{random.randint(100000, 999999)}"

def convert_fasta_to_csv(fasta_path, csv_path):
    records = parse_fasta(fasta_path)
    
    # Define headers: accession number, name, date, season, then x1 to x566
    headers = ["accession number", "name", "date", "season"] + [f"x{i}" for i in range(1, 567)]
    
    # Use QUOTE_ALL to match your sample format where every field is enclosed in quotes
    with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
        writer.writerow(headers)
        
        for name, seq in records:
            year = extract_year_from_name(name)
            epi_num = generate_epi()
            date_str = generate_random_date(year)
            season = str(year)
            
            # Prepare sequence characters, ensuring it matches 566 columns
            seq_chars = list(seq.lower())
            if len(seq_chars) < 566:
                seq_chars.extend([''] * (566 - len(seq_chars)))
            elif len(seq_chars) > 566:
                seq_chars = seq_chars[:566]
            
            row = [epi_num, name, date_str, season] + seq_chars
            writer.writerow(row)
            
    print(f"Successfully converted {len(records)} sequences to {csv_path}")

if __name__ == "__main__":
    # Replace 'sequences.fasta' and 'output.csv' with your actual file paths
    convert_fasta_to_csv('H3N2.fasta', 'H3N2_HA_HK_sequence.csv')
