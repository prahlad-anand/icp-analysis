"""
Initializes a SQLite database and loads all rows from cell-count.csv into cell_data.db.

Run:
    python load_data.py
"""

import csv
import os
import sqlite3

DB_PATH = "cell_data.db"
CSV_PATH = "cell-count.csv"

SCHEMA = """
CREATE TABLE IF NOT EXISTS subjects (
    subject_id  TEXT PRIMARY KEY,
    project_id  TEXT NOT NULL,
    condition   TEXT,
    age         INTEGER,
    sex         TEXT,
    treatment   TEXT,
    response    TEXT
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT NOT NULL,
    sample_type               TEXT,
    time_from_treatment_start INTEGER,
    b_cell                    INTEGER NOT NULL,
    cd8_t_cell                INTEGER NOT NULL,
    cd4_t_cell                INTEGER NOT NULL,
    nk_cell                   INTEGER NOT NULL,
    monocyte                  INTEGER NOT NULL,
    FOREIGN KEY (subject_id) REFERENCES subjects (subject_id)
);
"""


def init_db(conn: sqlite3.Connection):
    """
    Create tables.
    """
    conn.executescript(SCHEMA)
    conn.commit()


def load_csv(conn: sqlite3.Connection, csv_path: str=CSV_PATH):
    """
    Load the csv into the database, reading row-wise.
    """
    cursor = conn.cursor()
    subjects_seen = set()

    with open(csv_path, newline="", encoding="utf-8") as file:
        rows = csv.DictReader(file)
        for row in rows:
            project_id = row["project"]
            subject_id = row["subject"]
            sample_id  = row["sample"]

            # subjects table
            if subject_id not in subjects_seen:
                response = row["response"] if row["response"] else None
                cursor.execute(
                    """
                    INSERT OR IGNORE INTO subjects
                    (subject_id, project_id, condition, age, sex, treatment, response)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        subject_id,
                        project_id,
                        row["condition"],
                        int(row["age"]),
                        row["sex"],
                        row["treatment"],
                        response,
                    ),
                )
                subjects_seen.add(subject_id)

            # samples table
            cursor.execute(
                """
                INSERT OR IGNORE INTO samples
                (sample_id, subject_id, sample_type, time_from_treatment_start, b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    sample_id,
                    subject_id,
                    row["sample_type"],
                    int(row["time_from_treatment_start"]),
                    int(row["b_cell"]),
                    int(row["cd8_t_cell"]),
                    int(row["cd4_t_cell"]),
                    int(row["nk_cell"]),
                    int(row["monocyte"]),
                ),
            )

    conn.commit()


def main():
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
    conn = sqlite3.connect(DB_PATH)

    try:
        init_db(conn)
        print(f"Loading data from {CSV_PATH}\n")
        load_csv(conn)
        print(f"Database saved to {DB_PATH}\n")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
