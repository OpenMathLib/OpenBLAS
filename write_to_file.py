import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write contents to file.")
    parser.add_argument("contents", help="Contents.")
    parser.add_argument("file_path", help="File path.")

    args = parser.parse_args()

    f = open(args.file_path, "a")
    f.write(args.contents)
    f.close()
