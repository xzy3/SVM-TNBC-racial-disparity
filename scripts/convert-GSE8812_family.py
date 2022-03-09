#!/usr/bin/env python
import gzip
import csv
from contextlib import ExitStack
from Bio import Geo
from toolz.dicttoolz import valmap

with ExitStack() as stk:
    in_fd = stk.enter_context(gzip.open('data/GSE58812_family.soft.gz', 'rt'))
    out_fd = stk.enter_context(gzip.open('data/GSE58812.csv.gz', 'wt'))

    record_iter = Geo.parse(in_fd)
    # first two records are pure metadata and don't contain things we need
    next(record_iter)
    next(record_iter)
    gene_metadata_record = next(record_iter)
    gene_map = {
        row[0] : row[10] or row[11] or row[9]
            for row in gene_metadata_record.table_rows[1:]
    }

    header_row = [
        'sample',
        'accession',
        'age-at-diag',
        'death',
        'over-all-survival',
        'event-free-survival',
    ]
    header_row.extend(gene_map.values())

    writer = csv.DictWriter(out_fd, header_row)
    writer.writeheader()
    for rec in record_iter:
        print(rec)
        row = {
            'sample' : rec.entity_attributes['Sample_title'],
            'accession' : rec.entity_attributes['Sample_geo_accession'],
            'age-at-diag' : float(rec.entity_attributes['Sample_characteristics_ch1'][1].split(':')[-1]),
            'death' : int(rec.entity_attributes['Sample_characteristics_ch1'][4].split(':')[-1]),
            'over-all-survival' : int(rec.entity_attributes['Sample_characteristics_ch1'][5].split(':')[-1]),
            'event-free-survival' : int(rec.entity_attributes['Sample_characteristics_ch1'][3].split(':')[-1]),
        }
        row.update({
            gene_map[key] : float(val)
                for key,val in rec.table_rows[1:] })
        writer.writerow(valmap(str, row))
