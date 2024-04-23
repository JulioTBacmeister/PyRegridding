# Python Script: update_config.py
import yaml
from datetime import datetime, timedelta

def read_config_yaml(file_path):
    with open(file_path, "r") as config_file:
        config = yaml.safe_load(config_file)
    return config

def write_config_yaml(file_path, config):
    with open(file_path, "w") as config_file:
        yaml.dump(config, config_file)

def read_config(file_path):
    config = {}
    with open(file_path, 'r') as file:
        for line in file:
            # Skip blank lines
            if line.strip() == "":
                continue

            key, value = line.strip().split('=')
            if key in ('year','month','day','hour','Resubmit'):
                config[key] = int(value)
            else:
                config[key] = value
    return config

def write_config(file_path, config):
    with open(file_path, 'w') as file:
        for key, value in config.items():
            file.write(f"{key}={value}\n")


def increment_day(config, NoLeapYear=False ):
    # Create a datetime object from the config dictionary
    current_date = datetime(config['year'], config['month'], config['day'])
    
    # Increment the day using timedelta
    next_date = current_date + timedelta(days=1)

    # Manually adjust if next_date is February 29 and NoLeapYear is True
    if (NoLeapYear==True) and (next_date.month == 2) and (next_date.day == 29):
        # Adjust to March 1
        next_date = datetime(next_date.year, 3, 1)

    # Update the config dictionary with the new date
    config['year'], config['month'], config['day'] = next_date.year, next_date.month, next_date.day
    return config

def increment_hours(config, nhours=1, NoLeapYear=False ):
    # Create a datetime object from the config dictionary
    current_date = datetime(config['year'], config['month'], config['day'], config['hour'])
    
    # Increment the day using timedelta
    next_date = current_date + timedelta(hours=nhours)

    # Manually adjust if next_date is February 29 and NoLeapYear is True
    if (NoLeapYear==True) and (next_date.month == 2) and (next_date.day == 29):
        # Adjust to March 1
        next_date = datetime(next_date.year, 3, 1)

    # Update the config dictionary with the new date
    config['year'], config['month'], config['day'] , config['hour'] = next_date.year, next_date.month, next_date.day, next_date.hour
    return config

def increment_month(config):
    # Increment the month and handle month/year change if needed
    # This is a simplistic implementation and does not handle all edge cases
    config['month'] += 1
    if config['month'] > 12:
        config['month'] = 1
        config['year'] += 1
    return config

def decrement_Resubmit(config):
    # Decrement the Resubmit counter
    config['Resubmit'] -= 1
    return config

def main():
    file_path = './config.txt'  # Specify the path to your config file
    config = read_config(file_path)
    config = increment_day(config)
    write_config(file_path, config)

if __name__ == "__main__":
    main()
