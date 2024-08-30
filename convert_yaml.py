import yaml

# Specify the path to your YAML file and output file
input_file_path = 'sFree.yaml'
output_file_path = 'extracted_data.dat'

# Open and read the YAML file
with open(input_file_path, 'r') as file:
    documents = list(yaml.safe_load_all(file))

# Initialize a list to store the extracted data
extracted_data = []

# Iterate through each document and extract the values
for document in documents:
    for timestep, rows in document['data'].items():
        # Append the timestep and the three rows as a tuple
        for i in range(0, len(rows), 3):
            # Ensure there are three complete rows for each timestep
            if i + 2 < len(rows):
                extracted_data.append((timestep, rows[i][0], rows[i+1][0], rows[i+2][0]))

# Open the output file for writing
with open(output_file_path, 'w') as output_file:
    output_file.write("# Step DeltaE exp DeltaA\n")
    # Write the extracted data with whitespace delimiters
    for entry in extracted_data:
        output_file.write(f"{entry[0]} {entry[1]:.6f} {entry[2]:.6f} {entry[3]:.6f}\n")

print(f"Data successfully written to {output_file_path}")
