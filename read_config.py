#!/usr/bin/env python3

import argparse
from pathlib import Path


def read_config_file(file_path: Path) -> dict:
    config_data = {}
    with file_path.open("r") as file:
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
    parser = argparse.ArgumentParser(
        description="Read a config.h file and write to a specified build directory."
    )
    parser.add_argument("--file1", type=Path, help="Path to the config.h file.")
    parser.add_argument("--build_dir", type=Path, help="Path to the build directory.")

    args = parser.parse_args()

    config_data = read_config_file(args.file1)
    fdat = []
    for key, value in config_data.items():
        fdat.append(f"{key}={value}")
    result = "\n".join(fdat)

    # Ensure the build directory exists
    args.build_dir.mkdir(parents=True, exist_ok=True)

    # Write to the specified file in the build directory
    output_file_path = args.build_dir / "config.kconf"
    with output_file_path.open("a") as f:
        f.write(result)
