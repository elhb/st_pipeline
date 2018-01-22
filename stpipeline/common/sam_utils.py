""" 
This module contains some functions and utilities for ST SAM/BAM files
"""

import pysam
        
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



def get_annotations(bam_file_name, return_queue=None, return_filename=False):
    """
    Function that extracts the values of the XF tag in a bam file (gene_id from the annotation step)
    and returns a dictionary with counts for all values present
    """
    import subprocess
    import multiprocessing
    import sys

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

    if return_filename: annotations = (annotations,bam_file_name)
    if return_queue:
        return_queue.put(annotations)
    else:
        return annotations