#!/usr/bin/env python3

import argparse

def read_config_file(file_path: str) -> dict:
    config_data = {}
    with open(file_path, "r") as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith("#define"):
                parts = line.split()
                key = parts[1]
                if len(parts) == 3:
                    value = parts[2]
                    if value.isdigit():
                        value = int(value)
                    elif value.startswith('"') and value.endswith('"'):
                        value = value.strip('"')
                    config_data[key] = value
                elif len(parts) == 2:
                    config_data[key] = 1
    return config_data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read a config.h file.")
    parser.add_argument("file1", help="Path to the config.h file.")

    args = parser.parse_args()

    config_data = read_config_file(args.file1)
    fdat = []
    for key, value in config_data.items():
        fdat.append(f'{key}={value}')
    result = '\n'.join(fdat)

    f = open("./build/config.kconf", "a")
    f.write(result)
    f.close()
