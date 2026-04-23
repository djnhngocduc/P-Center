import os
import csv

def write_csv_rows(path: str, rows: list, fieldnames: list, append: bool = False):
    os.makedirs(os.path.dirname(path), exist_ok=True)

    mode = "a" if append else "w"
    file_exists = os.path.exists(path)

    with open(path, mode, newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)

        if (not append) or (not file_exists) or os.path.getsize(path) == 0:
            w.writeheader()

        for r in rows:
            w.writerow(r)