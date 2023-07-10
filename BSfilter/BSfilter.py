#!/usr/bin/env python3

import argparse
import datetime
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-fq', type=str)  # path/to/fastq/file
    parser.add_argument('-mC', type=float, default=1.0)  # methyl cytosine frequency allowed per read
    parser.add_argument('--count', type=int, default=100000)  # methyl cytosine count allowed per read

    args = parser.parse_args()

    start_time = datetime.datetime.now()
    fq_name, fq_ext = os.path.splitext(args.fq)
    new_fq_name = fq_name + "_filtered" + fq_ext

    with open(args.fq, "r") as handle, open(new_fq_name, "w") as fout:
        fq = list(handle)  # make subscriptable
        records = int(len(fq) / 4)  # number of records in the fastq file
        index_pos = 0

        for record in range(records):
            header_pos = index_pos  # the line position in the fasta file that is the header sequence
            seq_pos = index_pos + 1
            comment_pos = index_pos + 2
            phred_pos = index_pos + 3

            """Calculating m5C frequency"""
            mC_count = fq[seq_pos].count("C") + fq[seq_pos].count("c")  # number of cytosines in sequence
            N_count = len(fq[seq_pos].rstrip())  # length of the sequence
            mC_persistence = float(mC_count / N_count)  # percent of sequence that is cytosine

            """define filter and write to file"""
            if mC_persistence <= args.mC and mC_count <= args.count:
                fout.write("".join(f'{fq[header_pos]}'))
                fout.write("".join(f'{fq[seq_pos]}'))
                fout.write("".join(f'{fq[comment_pos]}'))
                fout.write("".join(f'{fq[phred_pos]}'))

            index_pos += 4

    with open(new_fq_name, 'r') as nfout:
        read_count2 = len(list(nfout)) / 4
        print("%i records in %s" % (read_count2, new_fq_name))

    end_time = datetime.datetime.now()
    print("time:", end_time - start_time, "\n")