import argparse


def merge_files(file1_path: str, file2_path: str) -> str:
    # Open files in read mode
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        content1 = file1.read()
        content2 = file2.read()
        merged_content = content1 + "\n" + content2
        return merged_content


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge and print content from two text files.")
    parser.add_argument("file1", help="Path to the first text file.")
    parser.add_argument("file2", help="Path to the second text file.")

    args = parser.parse_args()

    output = merge_files(args.file1, args.file2)
    print(output)
