def clean_fasta_headers(input_fasta, output_fasta):
    with open(input_fasta, 'r', encoding='utf-8') as infile, open(output_fasta, 'w', encoding='utf-8') as outfile:
        for line in infile:
            if line.startswith('>'):
                pipe_count = 0
                cut_index = -1
                for i, char in enumerate(line):
                    if char == '|':
                        pipe_count += 1
                        if pipe_count == 3:
                            cut_index = i
                            break
                if cut_index != -1:
                    line = line[:cut_index] + '\n'
            outfile.write(line)

if __name__ == "__main__":
    clean_fasta_headers('H3N2_HA_HK_sequence.fas', 'H3N2_HA_HK_sequence1.fas')
