#!/bin/bash
DB="/path/to/_broker_db/ztf_alerts_stream.db"
OUT="/path/to/dump.sql"

# 1. Dump schema
sqlite3 "$DB" ".schema" > "$OUT"

# 2. Loop over tables (now stripping \r characters)
sqlite3 "$DB" "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%';" \
| tr -d '\r' \
| while IFS= read -r table; do

  # Skip empty lines if any
  [[ -z "$table" ]] && continue

  echo "" >> "$OUT"
  echo "-- Data for $table" >> "$OUT"

  # Use the -cmd flag to set the mode before running the query
  sqlite3 "$DB" ".mode insert $table" "SELECT * FROM \"$table\" LIMIT 100;" >> "$OUT"

done

echo "COMMIT;" >> "$OUT"