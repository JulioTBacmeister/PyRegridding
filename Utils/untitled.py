import re



def extract_grid_values_from_file(file_path):
    # Read the code from the file
    with open(file_path, 'r') as file:
        code = file.read()

    # Regular expression to match grid values in elif statements in the gridInfo function
    pattern = r"elif\s+\(grid\s*==\s*'([^']+)'\)|elif\s+\(grid\s*==\s*'([^']+)'\)|elif\s+\(\(grid\s*==\s*'([^']+)'\)\s+or\s+\(grid\s*==\s*'([^']+)'\)\)"
    matches = re.findall(pattern, code)

    # Extract the grid values from matches
    grid_values = set()
    for match in matches:
        grid_values.update([val for val in match if val])  # Add non-empty values

    return sorted(grid_values)

# Provide the path to the module file here
file_path = 'path/to/your/module.py'
grid_values = extract_grid_values_from_file(file_path)

print("Possible grid values in gridInfo function:")
for val in grid_values:
    print(val)
