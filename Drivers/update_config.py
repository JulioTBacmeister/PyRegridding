# Python Script: update_config.py

"""
def read_config(file_path):
    config = {}
    with open(file_path, 'r') as file:
        for line in file:
            key, value = line.strip().split('=')
            config[key] = int(value)
    return config
"""
def read_config(file_path):
    config = {}
    with open(file_path, 'r') as file:
        for line in file:
            # Skip blank lines
            if line.strip() == "":
                continue

            key, value = line.strip().split('=')
            config[key] = int(value)
    return config

def write_config(file_path, config):
    with open(file_path, 'w') as file:
        for key, value in config.items():
            file.write(f"{key}={value}\n")

def increment_day(config):
    # Increment the day and handle month/year change if needed
    # This is a simplistic implementation and does not handle all edge cases
    config['day'] += 1
    if config['day'] > 31:
        config['day'] = 1
        config['month'] += 1
        if config['month'] > 12:
            config['month'] = 1
            config['year'] += 1
    return config

def increment_month(config):
    # Increment the month and handle month/year change if needed
    # This is a simplistic implementation and does not handle all edge cases
    config['month'] += 1
    if config['month'] > 12:
        config['month'] = 1
        config['year'] += 1
    return config

def main():
    file_path = './config.txt'  # Specify the path to your config file
    config = read_config(file_path)
    config = increment_day(config)
    write_config(file_path, config)

if __name__ == "__main__":
    main()
