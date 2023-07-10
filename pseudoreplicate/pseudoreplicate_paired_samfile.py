#!/usr/bin/env python3
import os
import time
import resource
import argparse
import pysam as pys
import numpy as np


if __name__ == "__main__":
    t1 = time.clock()

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='path to sam file')
    parser.add_argument('-o', type=str, default = "./", help='output directory')
    parser.add_argument('-r', type=int, default = 3, help='number of desired pseudoreplicates')
    parser.add_argument('-s', type = int, default=1, help='set seed')
    parser.add_argument('-b', action = 'store_true', help='output files in bam format, sorted by coordinate')
    parser.add_argument('--sorted', action = 'store_true', default=False, help='use if alignemnt file is already name sorted')

    args = parser.parse_args()

    in_file = args.i
    pseudo_replicates = args.r
    output_dir = args.o

    in_file_name, in_file_ext = os.path.splitext(in_file.split("/")[-1])
    out_file_name = os.path.join(output_dir,in_file_name )

    # Check file extension
    if not args.b and in_file_ext == ".bam":
        print('Your file must have a .sam extension.')
        print('If you are using a bam file, use -b flag.')
        exit(0)

    # Sort alignment file by name
    if not args.sorted:
        alignfile = str(out_file_name) + "_nsorted_tmp.sam"
        pys.sort("-o", alignfile, "-n", in_file)
    else: # need to convert to sam file
        alignfile = in_file

    # Open samfile, separate header lines, and create a set of qnames
    header = []
    read_count = 0

    with open(alignfile, 'r') as f:
        file = list(f)
        for line in file:
            if line[0] == "@":
                header.append(line)
            else:
                read_count +=1
        del file

    print('Your alignment file is being split...')

    # Create a new file for each replicate, and add the header
    split_file_list = []

    for i in range(pseudo_replicates):
        i += 1  # correct for 0-based indexing
        if args.b:
            name = os.path.join(output_dir,in_file_name.split("/")[-1]) + "_split" + str(i) + ".bam"
        else:
            name = os.path.join(output_dir,in_file_name.split("/")[-1]) + "_split" + str(i) + ".sam"
        split_file_list.append(name)
        with open(name, "w") as samfile_split:
            for j in header:
                samfile_split.write(j)

    # Add alignments to new file
    with open(alignfile, 'r') as f, open(split_file_list[0], 'a') as f1:
        samfile = list(f)[len(header):]     # start reading file at first read

        # Write the first read into a split file
        f1.write(samfile[0])

        # Initialize parameters for the rest of the reads
        index = 1  # The first read is at index 0. Starting at
        int = 0    # The first file containing reads is called _split1.sam
        # Begin partitioning reads into files randomly, keeping pairs together
        np.random.seed(args.s)  # doesn't work
        for i in range(read_count-1):
            if samfile[index].partition("\t")[0] == samfile[index-1].partition("\t")[0]:
                with open(split_file_list[int], 'a') as f2:
                    f2.write(samfile[index])
            else:
                int = np.random.randint(pseudo_replicates)
                with open(split_file_list[int], 'a') as f3:
                    f3.write(samfile[index])
            index += 1

    #convert output to bam is args.b, and remove all temp files
    if args.b:
        print('Converting split sams to bams and sorting by coordinate...')
        for file in split_file_list:
            name, ext = os.path.splitext(file)
            split_out_name = str(name) + '.bam'
            pys.sort("-o", split_out_name, file)

    if not args.sorted:
        os.remove(alignfile)

    t2 = time.clock()
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print("time: " + str(round(t2-t1, 3)) + " sec")
    print("memory: " + str(mem) + " bytes")