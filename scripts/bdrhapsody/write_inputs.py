#!/usr/bin/env python3

import argparse
from pathlib import Path


def parse_kv(values):
    parsed = []
    for item in values:
        key, value = item.split("=", 1)
        parsed.append((key, value))
    return parsed


def write_file_block(handle, field_name, paths):
    handle.write(f"{field_name}:\n")
    for path in paths:
        handle.write(" - class: File\n")
        handle.write(f'   location: "{path}"\n')


def write_single_file_block(handle, field_name, path):
    handle.write(f"{field_name}:\n")
    handle.write("  class: File\n")
    handle.write(f'  location: "{path}"\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--read", action="append", default=[])
    parser.add_argument("--file-field", action="append", default=[])
    parser.add_argument("--file-list-field", action="append", default=[])
    parser.add_argument("--scalar", action="append", default=[])
    parser.add_argument("--scalar-list", action="append", default=[])
    args = parser.parse_args()

    reads = [str(Path(path).resolve()) for path in args.read]
    file_fields = [(key, str(Path(value).resolve())) for key, value in parse_kv(args.file_field)]
    file_list_fields = [(key, str(Path(value).resolve())) for key, value in parse_kv(args.file_list_field)]
    scalars = parse_kv(args.scalar)
    scalar_lists = parse_kv(args.scalar_list)

    grouped_file_lists = {}
    for key, value in file_list_fields:
        grouped_file_lists.setdefault(key, []).append(value)

    grouped_scalar_lists = {}
    for key, value in scalar_lists:
        grouped_scalar_lists.setdefault(key, []).append(value)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("#!/usr/bin/env cwl-runner\n\n")
        handle.write("cwl:tool: rhapsody\n\n")
        write_file_block(handle, "Reads", reads)
        handle.write("\n")

        for key, value in file_fields:
            write_single_file_block(handle, key, value)
            handle.write("\n")

        for key, values in grouped_file_lists.items():
            write_file_block(handle, key, values)
            handle.write("\n")

        for key, value in scalars:
            handle.write(f"{key}: {value}\n")

        for key, values in grouped_scalar_lists.items():
            rendered = ", ".join(values)
            handle.write(f"{key}: [{rendered}]\n")


if __name__ == "__main__":
    main()
