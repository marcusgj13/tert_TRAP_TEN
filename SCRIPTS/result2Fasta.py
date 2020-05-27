from sys import argv

def parse_result_file(filename):
    """This function will read through a .result file as output
    by PROMALS and convert it to a dictionary."""
    seq_dict = {}

    with open(filename, 'r') as file_in:
        for line in file_in:
            if line[0] == '_':
                org_name = line.split()[0][1:]
                seq = line.split()[1]
                if org_name in seq_dict:
                    all_seq = seq_dict[org_name]
                    seq_dict[org_name] = all_seq + seq
                else:
                    seq_dict[org_name] = seq

    return seq_dict

def write_fasta_file(sequence_list, out_name):
    """Outputs a fasta file containing PROMALS alignments in a
    single line."""
    with open(out_name,'w') as file_out:
        #file_out.write('65 201\n')
        for key in sequence_list:
            file_out.write('>' + key + '\n')
            file_out.write(sequence_list[key] + '\n')


file_2_parse = argv[1]
new_filename = file_2_parse.split('.')[0] + '.fasta'
result_dict = parse_result_file(file_2_parse)
write_fasta_file(result_dict, new_filename)

