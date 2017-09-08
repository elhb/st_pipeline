""" 
This module contains some specific functionalities for
ST fastq files, mainly quality filtering functions.
"""

from stpipeline.common.utils import safeOpenFile, fileOk, is_fifo
from stpipeline.common.adaptors import removeAdaptor
from stpipeline.common.stats import qa_stats
import logging 
from itertools import izip
from sqlitedict import SqliteDict
import os
import re

def coroutine(func):
    """ 
    Coroutine decorator, starts coroutines upon initialization.
    """
    def start(*args, **kwargs):
        cr = func(*args, **kwargs)
        cr.next()
        return cr
    return start

def readfq(fp): # this is a generator function
    """ 
    Heng Li's fasta/fastq reader function.
    # https://github.com/lh3/readfq/blob/master/readfq.py
    # Unlicensed. 
    Parses fastq records from a file using a generator approach.
    :param fp: opened file descriptor
    :returns an iterator over tuples (name,sequence,quality)
    """
    cdef str last = '' # this is a buffer keeping the last unprocessed line
    cdef str l
    cdef str name
    cdef list seqs
    cdef int leng
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        #name, seqs, last = last[1:].partition(" ")[0], [], None
        name, seqs, last = last[1:], [], ''
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = ''
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break

@coroutine
def writefq(fp):  # This is a coroutine
    """ 
    Fastq writing generator sink.
    Send a (header, sequence, quality) triple to the instance to write it to
    the specified file pointer.
    """
    cdef str fq_format = '@{header}\n{sequence}\n+\n{quality}\n'
    cdef tuple record
    cdef str read
    try:
        while True:
            record = yield
            read = fq_format.format(header=record[0], sequence=record[1], quality=record[2])
            fp.write(read)
    except GeneratorExit:
        return
    
def quality_trim_index(bases, qualities, cutoff, base=33):
    """
    Function snippet and modified from CutAdapt 
    https://github.com/marcelm/cutadapt/
    
    Copyright (c) 2010-2016 Marcel Martin <marcel.martin@scilifelab.se>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN C

    Find the position at which to trim a low-quality end from a nucleotide sequence.

    Qualities are assumed to be ASCII-encoded as chr(qual + base).

    The algorithm is the same as the one used by BWA within the function
    'bwa_trim_read':
    - Subtract the cutoff value from all qualities.
    - Compute partial sums from all indices to the end of the sequence.
    - Trim sequence at the index at which the sum is minimal.
    
    This variant works on NextSeq data.
    With Illumina NextSeq, bases are encoded with two colors. 'No color' (a
    dark cycle) usually means that a 'G' was sequenced, but that also occurs
    when sequencing falls off the end of the fragment. The read then contains
    a run of high-quality G bases in the end.
    This routine works as the one above, but counts qualities belonging to 'G'
    bases as being equal to cutoff - 1.
    """
    s = 0
    max_qual = 0
    max_i = len(qualities)
    for i in reversed(xrange(max_i)):
        q = ord(qualities[i]) - base
        if bases[i] == 'G':
            q = cutoff - 1
        s += cutoff - q
        if s < 0:
            break
        if s > max_qual:
            max_qual = s
            max_i = i
    return max_i

def trim_quality(sequence,
                 quality,
                 min_qual=20, 
                 min_length=30, 
                 phred=33):    
    """
    Quality trims a fastq read using a BWA approach.
    It returns the trimmed record or None if the number of bases
    after trimming is below a minimum.
    :param sequence: the sequence of bases of the read
    :param quality: the quality scores of the read
    :param min_qual the quality threshold to trim (consider a base of bad quality)
    :param min_length: the minimum length of a valid read after trimming
    :param phred: the format of the quality string (33 or 64)
    :type sequence: str
    :type quality: str
    :type min_qual: integer
    :type min_length: integer
    :type phred: integer
    :return: A tuple (base, qualities) or (None,None)
    """
    if len(sequence) < min_length:
        return None, None
    # Get the position at which to trim (number of bases to trim)
    cut_index = quality_trim_index(sequence, quality, min_qual, phred)
    # Check if the trimmed sequence would have min length (at least)
    # if so return the trimmed read otherwise return None
    if (cut_index + 1) >= min_length:
        new_seq = sequence[:cut_index]
        new_qual = quality[:cut_index]
        return new_seq, new_qual
    else:
        return None, None
  
def check_umi_template(umi, template):
    """
    Checks that the UMI (molecular barcode) given as input complies
    with the pattern given in template.
    Returns True if the UMI complies
    :param umi: a molecular barcode
    :param template: a reg-based template with the same
                    distance of the UMI that should tell how the UMI should be formed
    :type umi: str
    :type template: str
    :return: True if the given molecular barcode fits the pattern given
    """
    p = re.compile(template)
    return p.match(umi) is not None


#TODO this approach uses too much memory
#     find a better solution (maybe Cython or C++)
def hashDemultiplexedReads(reads,
                           umi_start,
                           umi_end,
                           low_memory):
    """
    This function extracts the read name, the UMI and the x,y coordinates
    from the reads given as input (output of TaggD) and returns a dictionary
    with the clean read name as key and (x,y,umi) as value. X and Y correspond
    to the array coordinates of the barcode of the read and UMI is extracted from the read
    sequence.
    :param reads: path to a file with the fastq reads after demultiplexing (TaggD)
    :param umi_start: the start position of the UMI
    :param umi_end: the end position of the UMI
    :param low_memory: True to use a key-value db instead of dict
    :type reads: str
    :type umi_start: integer
    :type umi_end: integer
    :type low_memory: boolean
    :return: a dictionary of read_name -> (x,y,umi) tags where umi is optional
    """
    logger = logging.getLogger("STPipeline")
    
    if not os.path.isfile(reads):
        error = "Error parsing TaggD output, input file not present {}\n".format(reads)
        logger.error(error)
        raise RuntimeError(error)
    
    if low_memory:
        hash_reads = SqliteDict(autocommit=False, flag='c', journal_mode='OFF')
    else:
        hash_reads = dict()
    
    fastq_file = safeOpenFile(reads, "rU")
    for name, sequence, _ in readfq(fastq_file):
        # Assumes the header is like this
        # @NS500688:111:HNYW7BGXX:1:11101:13291:1099 1:N:0:TGCCCA B0:Z:GTCCCACTGGAACGACTGTCCCGCATC B1:Z:678 B2:Z:678
        header_tokens = name.split()
        assert(len(header_tokens) > 3)
        assert(len(sequence) >= umi_end)
        # Get the X and Y tags from the header of the read
        x = header_tokens[-2]
        y = header_tokens[-1]
        # The UMI is retrieved from the sequence
        umi = sequence[umi_start:umi_end]
        # We keep the same naming convention for the UMI attribute
        tags = (x,y,"B3:Z:{}".format(umi))
        # In order to save memory we truncate the read
        # name to only keep the unique part (lane, tile, x_pos, y_pos)
        # TODO this procedure is specific to only Illumina technology
        key = "".join(header_tokens[0].split(":")[-4:])
        hash_reads[key] = tags

    if low_memory: hash_reads.commit()
    fastq_file.close()    
    return hash_reads