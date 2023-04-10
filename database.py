import os
import sqlite3
import numpy as np
import io


def adapt_array(arr):
    """
    http://stackoverflow.com/a/31312102/190597 (SoulNibbler)
    """
    out = io.BytesIO()
    np.save(out, arr)
    out.seek(0)
    return sqlite3.Binary(out.read())


def convert_array(text):
    out = io.BytesIO(text)
    out.seek(0)
    return np.load(out)


# Converts np.array to TEXT when inserting
sqlite3.register_adapter(np.ndarray, adapt_array)

# Converts TEXT to np.array when selecting
sqlite3.register_converter("array", convert_array)


def create_db(path: str):
    con = sqlite3.connect(path, check_same_thread=False, timeout=100)
    cur = con.cursor()
    cur.execute("CREATE TABLE entropy(idx TEXT PRIMARY KEY, entropy REAL, sequence array)")
    con.commit()
    return con


class api:
    def __init__(self, path: str):
        if not os.path.isfile(path):
            con = create_db(path)
        else:
            con = sqlite3.connect(path, check_same_thread=False, timeout=100)
        self.con = con

    def insert(self, index: str, ent: float, seq: bytes):
        cur = self.con.cursor()
        cur.execute(
            "INSERT OR REPLACE INTO entropy VALUES (?, ?, ?)",
            (index, ent, seq)
        )
        self.con.commit()

    def bulk_insert(self, inserts: list):
        cur = self.con.cursor()
        cur.executemany(
            "INSERT OR REPLACE INTO entropy VALUES(?, ?, ?)",
            inserts
        )
        self.con.commit()

    def delete(self, index: str):
        cur = self.con.cursor()
        cur.execute(
            "DELETE FROM entropy WHERE idx = {}".format(
                index
            )
        )
        self.con.commit()

    def view(self):
        cur = self.con.cursor()
        res1 = cur.execute("SELECT idx, entropy FROM entropy WHERE entropy > 0")
        cur = self.con.cursor()
        res2 = cur.execute("SELECT idx, entropy FROM entropy WHERE entropy = 0")

        print("Number of entries with entropy 0: {}".format(len(list(res2))))
        print("Number of entries with entropy > 0: {}".format(len(list(res1))))
