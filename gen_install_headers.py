#!/usr/bin/env python3
import argparse
from pathlib import Path

def write_openblas_config_header(dest_dir, version, config_last_path, template_path):
    config_h_path = dest_dir / "openblas_config.h"
    with config_h_path.open('w') as f:
        f.write("#ifndef OPENBLAS_CONFIG_H\n")
        f.write("#define OPENBLAS_CONFIG_H\n")

        with config_last_path.open('r') as config_last:
            for line in config_last:
                if line.strip():
                    defines = line.split('#define ')
                    for define in defines:
                        if define.strip():
                            parts = define.split(maxsplit=1)
                            if len(parts) > 0:
                                macro_name = parts[0]
                                rest_of_line = " ".join(parts[1:]) if len(parts) > 1 else ""
                                line_to_write = f"#define OPENBLAS_{macro_name} {rest_of_line}"
                                f.write(f"{line_to_write.strip()}\n")
       

        f.write(f'#define OPENBLAS_VERSION " OpenBLAS {version} "\n')

        with template_path.open('r') as template:
            f.write(template.read())

        f.write("#endif /* OPENBLAS_CONFIG_H */\n")
    print(f"Generated openblas_config.h in {dest_dir}")

def write_f77blas_header(dest_dir, common_interface_path):
    f77blas_h_path = dest_dir / "f77blas.h"
    with f77blas_h_path.open('w') as f:
        f.write("#ifndef OPENBLAS_F77BLAS_H\n")
        f.write("#define OPENBLAS_F77BLAS_H\n")
        f.write('#include "openblas_config.h"\n')

        with common_interface_path.open('r') as common_interface:
            f.write(common_interface.read())

        f.write("#endif\n")
    print(f"Generated f77blas.h in {dest_dir}")

def write_cblas_header(dest_dir, cblas_path, symbol_prefix, symbol_suffix):
    cblas_h_path = dest_dir / "cblas.h"

    with cblas_path.open('r') as cblas_file:
        content = cblas_file.read()

    if symbol_prefix:
        content = re.sub(r'\bcblas', f'{symbol_prefix}cblas', content)
        content = re.sub(r'\bopenblas', f'{symbol_prefix}openblas', content)
        content = re.sub(f'{symbol_prefix}openblas_complex', 'openblas_complex', content)
        content = re.sub(r'\bgoto', f'{symbol_prefix}goto', content)

    if symbol_suffix:
        content = re.sub(r'\bcblas(\w*)', r'cblas\1' + symbol_suffix, content)
        content = re.sub(r'\bopenblas(\w*)', r'openblas\1' + symbol_suffix, content)
        content = re.sub(r'\bgoto(\w*)', r'goto\1' + symbol_suffix, content)
        content = re.sub(r'openblas_complex_(\w*)' + symbol_suffix, r'openblas_complex_\1', content)

    content = content.replace('common', 'openblas_config')

    with cblas_h_path.open('w') as f:
        f.write(content)

    print(f"Generated cblas.h in {dest_dir}")

def main():
    parser = argparse.ArgumentParser(description="Generate OpenBLAS headers")
    parser.add_argument('--dest-dir', required=True, help="Destination directory for headers")
    parser.add_argument('--version', required=True, help="OpenBLAS version")
    parser.add_argument('--config-last', required=True, help="Path to config_last.h")
    parser.add_argument('--template', required=True, help="Path to openblas_config_template.h")
    parser.add_argument('--common-interface', required=True, help="Path to common_interface.h")
    parser.add_argument('--cblas', required=True, help="Path to cblas.h")
    parser.add_argument('--symbol-prefix', default="", help="Symbol prefix for cblas.h")
    parser.add_argument('--symbol-suffix', default="", help="Symbol suffix for cblas.h")
    parser.add_argument('--generate-f77blas', action='store_true', help="Generate f77blas.h")
    parser.add_argument('--generate-cblas', action='store_true', help="Generate cblas.h")

    args = parser.parse_args()

    dest_dir = Path(args.dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    config_last_path = Path(args.config_last)
    template_path = Path(args.template)
    common_interface_path = Path(args.common_interface)
    cblas_path = Path(args.cblas)

    write_openblas_config_header(dest_dir, args.version, config_last_path, template_path)

    if args.generate_f77blas:
        write_f77blas_header(dest_dir, common_interface_path)

    if args.generate_cblas:
        write_cblas_header(dest_dir, cblas_path, args.symbol_prefix, args.symbol_suffix)

if __name__ == "__main__":
    main()
