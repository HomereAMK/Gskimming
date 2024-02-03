# Gskimming
 Analyses of Nuclear Reads Obtained using Genome skimming





## Pyscripts


# csv_to_tsv_stats.py
```bash
python csv_to_tsv_stats.py <input_csv_file> <output_tsv_file>
```

# badspecimens_identifier.py
```bash
# Current version : discarding individual with (error_rate) > 0.02) or (read_length != 'NA' and int(read_length) < 100)
python badspecimens_identifier.py <input.tsv> path/source/directory path/destination/directory
```


