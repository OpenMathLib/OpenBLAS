#!/usr/bin/env python3

import argparse
from pathlib import Path


def prepare_config_last(config_path, output_path, quad_precision, exprecision):
    config_path = Path(config_path)
    output_path = Path(output_path)

    with config_path.open("r") as config_file:
        config_lines = config_file.readlines()

    with output_path.open("w") as output_file:
        output_file.writelines(config_lines)
        if quad_precision:
            output_file.write("#define QUAD_PRECISION\n")
        if exprecision:
            output_file.write("#define EXPRECISION\n")


def main():
    parser = argparse.ArgumentParser(description="Prepare config_last.h from config.h")
    parser.add_argument("--config", required=True, help="Path to config.h")
    parser.add_argument("--output", required=True, help="Path to output config_last.h")
    parser.add_argument(
        "--quad-precision", action="store_true", help="Enable QUAD_PRECISION"
    )
    parser.add_argument("--exprecision", action="store_true", help="Enable EXPRECISION")

    args = parser.parse_args()

    prepare_config_last(args.config, args.output, args.quad_precision, args.exprecision)


if __name__ == "__main__":
    main()
