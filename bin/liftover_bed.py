#!/usr/bin/env python
"""
liftover_bed.py: BED4 liftover using the Python liftover library.
Argument order: input.bed chain chrom_sizes output.bed unmapped.bed
chrom_sizes is a .fai (or plain chrom<TAB>size) file for the target genome;
lifted intervals falling outside a chromosome's real length are rejected.
If chrom_sizes is missing, empty, or unparseable (e.g. the NO_FILE_FASTA_FAI
placeholder used when no .fai was supplied), the bounds check is skipped
and only the pre-existing chrom/strand consistency check applies.
"""
import os
import sys
from liftover import get_lifter


def read_chrom_sizes(path):
    sizes = {}
    if not path or not os.path.isfile(path):
        return sizes
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 2:
                continue
            try:
                sizes[cols[0]] = int(cols[1])
            except ValueError:
                continue
    return sizes


def main():
    if len(sys.argv) != 6:
        sys.exit("Usage: liftover_bed.py <input.bed> <chain.gz> <chrom_sizes> <output.bed> <unmapped.bed>")

    in_bed, chain_path, chrom_sizes_path, out_bed, unmapped_bed = sys.argv[1:6]

    converter = get_lifter(chain_path)
    chrom_sizes = read_chrom_sizes(chrom_sizes_path)
    lifted = []
    unmapped = []

    with open(in_bed) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 4:
                unmapped.append(line)
                continue

            chrom, start, end, name = cols[0], int(cols[1]), int(cols[2]), cols[3]

            try:
                start_hits = converter[chrom][start]
                end_hits   = converter[chrom][end - 1]  # BED end is exclusive
            except Exception:
                unmapped.append(line)
                continue

            if start_hits and end_hits:
                s_chrom, s_pos, s_strand = start_hits[0]
                e_chrom, e_pos, e_strand = end_hits[0]
                if s_chrom == e_chrom and s_strand == e_strand:
                    if s_strand == '+':
                        new_start, new_end = s_pos, e_pos + 1
                    else:
                        new_start, new_end = e_pos, s_pos + 1

                    chrom_len = chrom_sizes.get(s_chrom)
                    out_of_bounds = new_start < 0 or new_start >= new_end
                    if chrom_sizes:
                        # only enforce the length check when sizes were actually available,
                        # an unrecognised chromosome is itself a reason to reject
                        out_of_bounds = out_of_bounds or chrom_len is None or new_end > chrom_len
                    if out_of_bounds:
                        unmapped.append(line)
                    else:
                        lifted.append('\t'.join([s_chrom, str(new_start), str(new_end), name]))
                else:
                    unmapped.append(line)
            else:
                unmapped.append(line)

    with open(out_bed, 'w') as fh:
        for r in lifted:
            fh.write(r + '\n')

    with open(unmapped_bed, 'w') as fh:
        for line in unmapped:
            fh.write(line + '\n')


if __name__ == '__main__':
    main()
