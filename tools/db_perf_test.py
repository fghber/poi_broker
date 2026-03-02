import sqlite3
import os
import time

def main():
    # Path to the database
    db_path = r"c:\dev\python\_broker_db\ztf_alerts_stream.db"

    # File existence and size
    exists = os.path.exists(db_path)
    print("exists:", exists)
    print("size_mb:", round(os.path.getsize(db_path) / 1024 / 1024, 2) if exists else None)

    # Connect to the database
    con = sqlite3.connect(db_path)
    cur = con.cursor()

    # List tables
    cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = [r[0] for r in cur.fetchall()]
    print("tables:", tables[:10])

    # Index list
    cur.execute("PRAGMA index_list('featuretable')")
    idx = cur.fetchall()
    print("index_count:", len(idx))
    print("indexes:", idx[:8])

    # Count rows
    t = time.perf_counter()
    cur.execute("SELECT COUNT(*) FROM featuretable")
    n = cur.fetchone()[0]
    t1 = time.perf_counter()
    print("rows:", n, "count_s:", round(t1 - t, 3))

    # Top 100 query timing
    t = time.perf_counter()
    cur.execute("SELECT date_alert_mjd, alert_id FROM featuretable ORDER BY date_alert_mjd DESC LIMIT 100")
    cur.fetchall()
    t1 = time.perf_counter()
    print("top100_s:", round(t1 - t, 3))

    # Filtered count timing
    t = time.perf_counter()
    cur.execute("SELECT COUNT(*) FROM featuretable WHERE ztf_object_id IS NOT NULL")
    cur.fetchone()
    t1 = time.perf_counter()
    print("filtered_count_s:", round(t1 - t, 3))

    # Close connection
    con.close()

if __name__ == "__main__":
    main()

"""
Cold Start (after DB creation, no indexes):
$>python db_perf_test.py
exists: True
size_mb: 3058.75
tables: ['featuretable', 'crossmatches', 'classification']
index_count: 1
indexes: [(0, 'sqlite_autoindex_featuretable_1', 1, 'pk', 0)]
rows: 1363178 count_s: 1.633
top100_s: 0.409
filtered_count_s: 33.522
"""

"""
Without Indexes:
$>python db_perf_test.py
exists: True
size_mb: 3058.75
tables: ['featuretable', 'crossmatches', 'classification']
index_count: 1
indexes: [(0, 'sqlite_autoindex_featuretable_1', 1, 'pk', 0)]
rows: 1363178 count_s: 0.062
top100_s: 0.453
filtered_count_s: 2.044
"""

"""
With Indexes:
python db_perf_test.py
exists: True
size_mb: 3303.07
tables: ['featuretable', 'crossmatches', 'classification', 'sqlite_stat1', 'sqlite_stat4']
index_count: 8
indexes: [(0, 'idx_featuretable_locus_dec', 0, 'c', 0), (1, 'idx_featuretable_locus_ra', 0, 'c', 0), (2, 'idx_featuretable_ant_passband', 0, 'c', 0), (3, 'idx_featuretable_locus_id', 0, 'c', 0), (4, 'idx_featuretable_ztf_object_id', 0, 'c', 0), (5, 'idx_featuretable_alert_id', 0, 'c', 0), (6, 'idx_featuretable_date_alert_mjd', 0, 'c', 0), (7, 'sqlite_autoindex_featuretable_1', 1, 'pk', 0)]
rows: 1363178 count_s: 0.016
top100_s: 0.001
filtered_count_s: 0.032
"""