#!/usr/bin/env python
import lzma
import csv
import sys
from contextlib import ExitStack
from toolz.dicttoolz import valmap,valfilter
from Bio import Geo

with ExitStack() as stk:
    in_fd = stk.enter_context(lzma.open('data/GSE142102_family.soft.xz', 'rt'))
    out_fd = stk.enter_context(lzma.open('data/GSE142102.csv.xz', 'wt', preset=9))

    import ipdb
    ipdb.set_trace()
    record_iter = Geo.parse(in_fd)
    # first two records are pure metadata and don't contain things we need
    next(record_iter)
    next(record_iter)
    gene_metadata_record = next(record_iter)
    gene_map = {
        row[0] : row[7].split('//', 2)[1].strip()
            for row in gene_metadata_record.table_rows[1:] if row[7] != '---'
    }

    header_row = [
        'sample',
        'accession',
        'age-at-diag',
        'death',
        'batch'
    ]
    header_row.extend(gene_map.values())

    writer = csv.DictWriter(out_fd, header_row)
    writer.writeheader()
    for rec in record_iter:
        print(rec)
        vital_status = None
        age = None
        batch = None
        for val in rec.entity_attributes['Sample_characteristics_ch1']:
            try:
                if val.startswith('vital status:'):
                    vital_status = val

                elif val.startswith('batch:'):
                    batch = int(val.split(':')[1])

                elif val.startswith('age:'):
                    age = int(val.split(':')[1])

            except ValueError:
                # drop NA values
                pass

        if vital_status is None or age is None or batch is None:
            continue

        if vital_status != 'vital status: Dead' and vital_status != 'vital status: Alive':
            print('Unknown vital status', vital_status)

        row = {
            'sample' : rec.entity_attributes['Sample_title'],
            'accession' : rec.entity_attributes['Sample_geo_accession'],
            'death' : vital_status == 'vital status: Dead',
            'age-at-diag' : age,
            'batch' : batch
        }
        row.update({
            gene_map[key] : float(val)
                for key,val in rec.table_rows[1:] if key in gene_map
        })
        try:
            writer.writerow(valmap(str, row))
        except ValueError:
            ipdb.set_trace()
            pass

