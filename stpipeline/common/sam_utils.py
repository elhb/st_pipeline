"""
This module contains some functions and utilities for ST SAM/BAM files
"""

import pysam
import subprocess
import multiprocessing
import sys

def convert_to_AlignedSegment(header, sequence, quality,
                              barcode_sequence, umi_sequence):
    """
    This function converts the input variables
    (header,sequence,quality,barcode_sequence,umi_sequence)
    to a unaligned pysam.AlignedSegment with the umi and barcode
    informations as the following tags:
        Tag  Value
        "B0" barcode_sequence
        "B3" umi_sequence
    :param header: string with the header information
    :param sequence: string with the DNA/RNA sequence
    :param quality: string with the base calling quality values
    :param barcode_sequence: string with the barcode sequence
    :param umi_sequence: string with the unique molecular identifier sequence
    """

    # create
    aligned_segment = pysam.AlignedSegment()

    # Set the standard values
    # Header must not contain empty spaces
    aligned_segment.query_name = header.split()[0]
    aligned_segment.query_sequence = sequence
    aligned_segment.query_qualities = pysam.qualitystring_to_array(quality)

    # setting the flag to un_mapped
    aligned_segment.flag |= pysam.FUNMAP

    # Set the tags
    aligned_segment.set_tag('B0', barcode_sequence)
    aligned_segment.set_tag('B3', umi_sequence)
    aligned_segment.set_tag('RG', '0')

    return aligned_segment

def merge_bam(merged_file_name, files_to_merge, ubam=False, samtools=True):
    """
    Function for merging partial bam files after annotation.
    also counts the number of reads for different types of annotations (the XF tags of the reads)
    and returns these counts as a dict: anotation=>count
    :param merged_file_name: name of the merged output bamfile
    :param files_to_merge: list with names of the partial bamfiles to merge
    :param ubam: indicates unaligned bam file (True or False, default False)
    :returns: the number of annotated records
    """
    annotations = {}

    if samtools: # use samtools to merge the bamfile instead of pysam, should be faster

        if ubam: # the samtools merge version is not implemented for unaligned bam files
            sys.stderr.write('The samtools merge version is not implemented for unaligned bam files.\n')
            raise NotImplementedError

        # Start to merge the files
        threads=len(files_to_merge)
        arguments = ['samtools','merge','-@',str(threads), merged_file_name]
        for filename in files_to_merge: arguments.append(filename)
        samtools_merge = subprocess.Popen(
            arguments,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

        # get the annotations in parallel
        keyword_arguments = { 'return_queue':multiprocessing.Queue() }
        sub_processes = [
                multiprocessing.Process(target=get_annotations, args=[bam_file_name], kwargs=keyword_arguments)
                for bam_file_name in files_to_merge
            ]
        for process in sub_processes: process.start()
        for process in sub_processes: process.join()
        while not keyword_arguments['return_queue'].empty():
            temp_dict = keyword_arguments['return_queue'].get()
            for annotation,count in temp_dict.iteritems():
                try:
                    annotations[annotation] += count
                except KeyError:
                    annotations[annotation] = count

        # check the completion of the samtools merge sub process
        stdout, stderr = samtools_merge.communicate()
        if len(stderr) > 0 or samtools_merge.returncode != 0:
            msg = "The samtools merge subprocess generated an error while running the stpipeline.common.sam_utils.merge_bam function (exit code {}).\n{}\n".format(samtools_merge.returncode,stderr)
            sys.stderr.write(msg)
            sys.exit(1)

    else:
        with pysam.AlignmentFile(files_to_merge[0], mode='rb', check_sq=(not ubam)) as input_bamfile:
            merged_file = pysam.AlignmentFile(merged_file_name,
                                              mode="wb", template=input_bamfile)
        # Simply merges the BAM files and creates a counter of annotated records
        for file_name in files_to_merge:
            input_bamfile = pysam.AlignmentFile( file_name, mode='rb', check_sq=(not ubam) )
            for record in input_bamfile.fetch(until_eof=True):
                merged_file.write(record)
                if ubam: annotation = None
                else: annotation = record.get_tag("XF")
                try:
                    annotations[annotation] += 1
                except KeyError:
                    annotations[annotation] = 1

    return sum(annotations.values())

def get_annotations(bam_file_name, return_queue=None):
    """
    Function that extracts the values of the XF tag in a bam file (gene_id from the annotation step)
    and returns a dictionary with counts for all values present
    """
    annotations = dict()

    samtools_view = subprocess.Popen(
        ['samtools','view',bam_file_name],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)

    grep = subprocess.Popen(
        ['grep','-E','XF\:Z\:.+(\t|$)','--only-matching'],
        stdin=samtools_view.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)

    unique = subprocess.Popen(
        ['uniq','-c'],
        stdin=grep.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)

    while True:
        line = unique.stdout.readline().rstrip().lstrip().split()
        if line:
            try:
                annotations[line[1].split('XF:Z:')[1]] += int(line[0])
            except KeyError:
                annotations[line[1].split('XF:Z:')[1]] = int(line[0])
        else: break

    errmsg = "A {} subprocess generated an error while running the stpipeline.common.sam_utils.get_annotations function.\n{}\n"
    stdout, stderr = samtools_view.communicate()
    if len(stderr) > 0:
        sys.stderr.write(msg.format('samtools view',stderr))
        sys.exit(1)

    stdout, stderr = grep.communicate()
    if len(stderr) > 0:
        sys.stderr.write(msg.format('grep',stderr))
        sys.exit(1)

    stdout, stderr = unique.communicate()
    if len(stderr) > 0:
        sys.stderr.write(msg.format('uniq',stderr))
        sys.exit(1)

    if return_queue:
        return_queue.put(annotations)
    else:
        return annotations
