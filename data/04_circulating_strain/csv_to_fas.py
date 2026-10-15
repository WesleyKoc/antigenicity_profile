import csv

def csv_to_fasta(csv_path, fasta_path):
    with open(csv_path, 'r', encoding='utf-8') as csvfile, open(fasta_path, 'w', encoding='utf-8') as fastafile:
        reader = csv.reader(csvfile)
        
        # Read the header row to skip it
        header_row = next(reader, None)
        if not header_row:
            print("The CSV file is empty.")
            return
        
        count = 0
        for row in reader:
            if not row or len(row) < 4:
                continue
                
            accession = row[0]
            name = row[1]
            date = row[2]
            
            # Replace spaces in the name with underscores to match your target style
            formatted_name = name.replace(' ', '_')
            
            # Extract sequence characters from columns x1 onwards (index 4 to the end)
            # Filter out any empty strings if padding was applied
            seq_chars = [char for char in row[4:] if char != '']
            sequence = "".join(seq_chars).upper()
            
            # Construct the header in the requested format: >accession_|_date_|_name
            fasta_header = f">{accession}_|_{date}_|_{formatted_name}.upper()"
            
            # Write to FASTA file
            fastafile.write(f"{fasta_header}\n{sequence}\n")
            count += 1
            
    print(f"Successfully converted {count} rows from {csv_path} to {fasta_path}")

if __name__ == "__main__":
    # Replace 'output.csv' and 'sequences.fas' with your actual file names/paths
    csv_to_fasta('H3N2_HA_HK_sequence.csv', 'H3N2_HA_HK_sequence.fas')
