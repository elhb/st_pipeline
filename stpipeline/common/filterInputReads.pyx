import multiprocessing
import sys
import time
from stpipeline.common.utils import safeOpenFile, fileOk, is_fifo
from stpipeline.common.fastq_utils import *
from stpipeline.common.adaptors import removeAdaptor
from stpipeline.common.stats import qa_stats
import logging 
from itertools import izip
from sqlitedict import SqliteDict
import os
import re
import Queue

class InputReadsFilter():
    """
    DESC
    """
    
    def __init__(self, verbose=True, max_reads_in_memory=10000000, stat_line_interwall=100000, chunk_size=100000):
        self.verbose = verbose
        self.aborted = False
        self.stat_line_interwall = stat_line_interwall
        self.chunk_size = chunk_size
        self.max_chunks_in_queue = min(32767,max_reads_in_memory/self.chunk_size*2)
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: initiation completed.\n')
    
    def input_arguments(
                    self,
                    fw, 
                    rv,
                    out_fw,
                    out_rw,
                    out_rw_discarded,
                    filter_AT_content,
                    filter_GC_content,
                    umi_start, 
                    umi_end,
                    min_qual, 
                    min_length,
                    polyA_min_distance, 
                    polyT_min_distance, 
                    polyG_min_distance, 
                    polyC_min_distance,
                    polyN_min_distance,
                    qual64,
                    umi_filter,
                    umi_filter_template,
                    umi_quality_bases,
                    adaptor_missmatches):

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: Getting input arguments.\n')

        self.logger = logging.getLogger("STPipeline")
        if not (os.path.isfile(fw) or is_fifo(fw)) or not (os.path.isfile(rv) or is_fifo(rv)):
            error = "Error doing quality trimming, input file/s not present {}\n{}\n".format(fw,rv)
            logger.error(error)
            raise RuntimeError(error)
        
        # Check if discarded files must be written out 
        self.keep_discarded_files = out_rw_discarded is not None
        
        # Build fake sequence adaptors with the parameters given
        self.adaptorA = "".join("A" for k in xrange(polyA_min_distance))
        self.adaptorT = "".join("T" for k in xrange(polyT_min_distance))
        self.adaptorG = "".join("G" for k in xrange(polyG_min_distance))
        self.adaptorC = "".join("C" for k in xrange(polyC_min_distance))
        self.adaptorN = "".join("N" for k in xrange(polyN_min_distance))
        
        # Not recommended to do adaptor trimming for adaptors smaller than 5
        self.do_adaptorA = polyA_min_distance >= 5
        self.do_adaptorT = polyT_min_distance >= 5
        self.do_adaptorG = polyG_min_distance >= 5
        self.do_adaptorC = polyC_min_distance >= 5
        self.do_adaptorN = polyN_min_distance >= 5
        self.do_AT_filter = filter_AT_content > 0
        self.do_GC_filter = filter_GC_content > 0
        
        # Quality format
        self.phred = 64 if qual64 else 33

        self.umi_filter = umi_filter
        self.umi_start = umi_start
        self.umi_end = umi_end
        self.umi_filter_template = umi_filter_template
        self.umi_quality_bases = umi_quality_bases
        self.filter_AT_content = filter_AT_content
        self.filter_GC_content = filter_GC_content
        self.adaptor_missmatches = adaptor_missmatches
        self.min_length = min_length
        self.min_qual = min_qual

        self.fw = fw
        self.rv = rv
        self.out_fw = out_fw
        self.out_rw = out_rw
        self.out_rw_discarded = out_rw_discarded
    
    def run(self, ):
        """
        Function that starts the subprocesses reading the file(s)
        """
        
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process pid='+str(os.getpid())+'.\n')
        
        # create queues and connections
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - Creating queues and connections.\n')
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - Creating input reads queue.\n')
        self.input_read_queue = multiprocessing.Queue()
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - Creating output reads queue.\n')
        self.output_read_queue = multiprocessing.Queue()
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - Creating counter pipe.\n')
        self.counter_connection_send_end, counter_connection_recv_end = multiprocessing.Pipe()
        
        import ctypes
        self.reader_running = multiprocessing.Value(ctypes.c_bool, True)
        self.workers_running= multiprocessing.Value(ctypes.c_bool, True)
        
        # start reader subprocess
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - Starting reader.\n')
        self.reader = multiprocessing.Process(target=self.input_files_reader)
        self.reader.start()

        # start worker pool
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - Starting pool.\n')
        self.worker_pool = [multiprocessing.Process(target=self.parallel_worker_function) for i in range(multiprocessing.cpu_count()-2)]
        for process in self.worker_pool: process.start()

        # start writer subprocess
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - Creating writer.\n')
        self.writer = multiprocessing.Process(target=self.output_files_writer)
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - Stariting writer.\n')
        self.writer.start()
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - writer started.\n')

        # wait for reader to complete
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - waiting for reader.\n')
        self.reader.join()
        
        # wait for worker pool to complete
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - waiting for pool.\n')
        for process in self.worker_pool: process.join()
        self.workers_running.value = False
        
        # wait for writer to complete
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - waiting for writer.\n')
        self.writer.join()
        
        # get the counters from the writer process
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - fetching counter.\n')
        read_pair_counters = counter_connection_recv_end.recv()
        counter_connection_recv_end.close()
        self.counter_connection_send_end.close()

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - updating logs and checking for prescence of files.\n')
        # Write info to the log
        self.logger.info("Trimming stats total reads (pair): {}".format(read_pair_counters['total_reads']))
        self.logger.info("Trimming stats {} reads have been dropped!".format(read_pair_counters['dropped_rv'])) 
        perc2 = '{percent:.2%}'.format(percent= float(read_pair_counters['dropped_rv']) / float(read_pair_counters['total_reads']) )
        self.logger.info("Trimming stats you just lost about {} of your data".format(perc2))
        self.logger.info("Trimming stats reads remaining: {}".format(read_pair_counters['total_reads'] - read_pair_counters['dropped_rv']))
        self.logger.info("Trimming stats dropped pairs due to incorrect UMI: {}".format(read_pair_counters['dropped_umi_template']))
        self.logger.info("Trimming stats dropped pairs due to low quality UMI: {}".format(read_pair_counters['dropped_umi']))
        self.logger.info("Trimming stats dropped pairs due to high AT content: {}".format(read_pair_counters['dropped_AT']))
        self.logger.info("Trimming stats dropped pairs due to high GC content: {}".format(read_pair_counters['dropped_GC']))
        self.logger.info("Trimming stats dropped pairs due to presence of artifacts: {}".format(read_pair_counters['dropped_adaptor']))
        
        # Check that output file was written ok
        if not fileOk(self.out_rw):
            error = "Error doing quality trimming checks of raw reads." \
            "\nOutput file not present {}\n".format(self.out_rw)
            logger.error(error)
            raise RuntimeError(error)
        
        # Adding stats to QA Stats object
        qa_stats.input_reads_forward = read_pair_counters['total_reads']
        qa_stats.input_reads_reverse = read_pair_counters['total_reads']
        qa_stats.reads_after_trimming_forward = (read_pair_counters['total_reads'] - read_pair_counters['dropped_rv'])
        qa_stats.reads_after_trimming_reverse = (read_pair_counters['total_reads'] - read_pair_counters['dropped_rv'])
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: main process - completed.\n')
        
    def input_files_reader(self, ):

        cdef double start_time = time.time()
        cdef int count = 0
        cdef double last_time = start_time
        cdef int last_count = count
        cdef str identity = 'READER'

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: Starting input files reader pid='+str(os.getpid())+'.\n')

        # Open fastq files with the fastq parser
        fw_file = safeOpenFile(self.fw, "rU")
        rv_file = safeOpenFile(self.rv, "rU")

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: input files reader - starting to parse reads.\n')
        if self.verbose: self.print_stat_line(start_time, last_time, count, last_count, identity, header=True)
        cdef list chunk = []
        cdef tuple read_one
        cdef tuple read_two
        for read_one, read_two in izip(readfq(fw_file), readfq(rv_file)):
            chunk.append((read_one, read_two))
            if len(chunk) == self.chunk_size:
                self.input_read_queue.put( chunk )
                chunk = []
            count += 1
            if self.verbose and count % self.stat_line_interwall == 0:
                self.print_stat_line(start_time, last_time, count, last_count, identity)
                last_time = time.time()
                last_count = count

        self.input_read_queue.put( chunk )
        
        if self.verbose:
            self.print_stat_line(start_time, last_time, count, last_count, identity)
            last_time = time.time()
            last_count = count

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: input files reader - all reads parsed.\n')
        fw_file.close()
        rv_file.close()
        
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: input files reader - closing and joining queue.\n')
        self.input_read_queue.close()
        self.input_read_queue.join_thread()
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: input files reader - shutting down.\n')
        self.reader_running.value = False

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: input files reader - completed.\n')

    def output_files_writer(self, ):

        cdef double start_time = time.time()
        cdef int count = 0
        cdef double last_time = start_time
        cdef int last_count = count
        cdef str identity = 'WRITER'

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: Starting output files writer pid='+str(os.getpid())+'.\n')
        
        # Create output file writers
        out_rv_handle = safeOpenFile(self.out_rw, 'w')
        out_rv_writer = writefq(out_rv_handle)
        out_fw_handle = safeOpenFile(self.out_fw, 'w')
        out_fw_writer = writefq(out_fw_handle)
        if self.keep_discarded_files:
            out_rv_handle_discarded = safeOpenFile(self.out_rw_discarded, 'w')
            out_rv_writer_discarded = writefq(out_rv_handle_discarded)
        
        # Some counters
        cdef dict read_pair_counters = {
            'total_reads':0,
            'dropped_rv':0,
            'dropped_umi':0,
            'dropped_umi_template':0,
            'dropped_AT': 0,
            'dropped_GC':0,
            'dropped_adaptor':0
        }
        cdef list chunk
        cdef dict read_pair

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: output files writer - starting to parse reads from queue.\n')
        while self.workers_running.value or not self.output_read_queue.empty():
        
            while not self.output_read_queue.empty():
                
                chunk = self.output_read_queue.get()
                for read_pair in chunk:
        
                    read_pair_counters['total_reads'] += 1
                    if read_pair['discard_reson'] == 'dropped_umi_template': read_pair_counters['dropped_umi_template'] += 1
                    if read_pair['discard_reson'] == 'dropped_umi':          read_pair_counters['dropped_umi'] += 1
                    if read_pair['discard_reson'] == 'dropped_AT':           read_pair_counters['dropped_AT'] += 1
                    if read_pair['discard_reson'] == 'dropped_GC':           read_pair_counters['dropped_GC'] += 1
                    if read_pair['discard_reson'] == 'dropped_adaptor':      read_pair_counters['dropped_adaptor'] += 1                
                    
                    # Write reverse read to output
                    if not read_pair['discard_reson']:
                        out_rv_writer.send((read_pair['header_rv'], read_pair['sequence_rv'], read_pair['quality_rv']))
                        out_fw_writer.send((read_pair['header_fw'], read_pair['sequence_fw'], read_pair['quality_fw']))
                    else:
                        read_pair_counters['dropped_rv'] += 1  
                        if self.keep_discarded_files:
                            out_rv_writer_discarded.send((read_pair['header_rv'], read_pair['orig_sequence_rv'], read_pair['orig_quality_rv']))
        
                    count += 1
                    if self.verbose and count % self.stat_line_interwall == 0:
                        self.print_stat_line(start_time, last_time, count, last_count, identity)
                        last_time = time.time()
                        last_count = count
            #if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: output files writer - queue empty.\n')

        if self.verbose:
            self.print_stat_line(start_time, last_time, count, last_count, identity)
            last_time = time.time()
            last_count = count

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: output files writer - all reads parsed.\n')

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: output files writer - sending counter to parent process.\n')
        self.counter_connection_send_end.send(read_pair_counters)
        
        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: output files writer - closing outfiles.\n')
        out_rv_handle.flush()
        out_rv_handle.close()
        out_rv_writer.close()
        out_fw_handle.flush()
        out_rv_writer.close()
        out_fw_writer.close()
        if self.keep_discarded_files:
            out_rv_handle_discarded.flush()
            out_rv_handle_discarded.close()
            out_rv_writer_discarded.close()

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: output files writer - completed.\n')

    def parallel_worker_function(self, ):

        start_time = time.time()
        count = 0
        last_time = start_time
        last_count = count
        identity = 'WORKER_'+str(os.getpid())

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: Starting parallel worker process pid='+str(os.getpid())+'.\n')

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: parallel worker process pid='+str(os.getpid())+' - starting to parse reads.\n')
        while self.reader_running.value or not self.input_read_queue.empty():
            
            while not self.input_read_queue.empty():
                
                try:
                    in_chunk = self.input_read_queue.get(timeout=10)
                except Queue.Empty: continue
                out_chunk = []
                for (header_fw, sequence_fw, quality_fw), (header_rv, sequence_rv, quality_rv) in in_chunk:
                    
                    read_pair = {
                        'header_fw':header_fw,
                        'sequence_fw':sequence_fw,
                        'quality_fw':quality_fw,
                        'header_rv':header_rv,
                        'sequence_rv':sequence_rv,
                        'quality_rv':quality_rv,
                        'orig_sequence_rv':sequence_rv,
                        'orig_quality_rv':quality_rv,
                        'discard_reson':None
                    }
                    
                    if not read_pair['sequence_fw'] or not read_pair['sequence_rv']:
                        error = "Error doing quality trimming, Checks of raw reads.\n" \
                        "The input files are not of the same length"
                        #"The input files {},{} are not of the same length".format(fw,rv)
                        logger.error(error)
                    
                    if read_pair['header_fw'].split()[0] != read_pair['header_rv'].split()[0]:
                        logger.warning("Pair reads found with different " \
                                       "names {} and {}".format(read_pair['header_fw'],read_pair['header_rv']))
                                
                    # If we want to check for UMI quality and the UMI is incorrect
                    # then we discard the reads
                    if self.umi_filter \
                    and not check_umi_template(read_pair['sequence_fw'][self.umi_start:self.umi_end], self.umi_filter_template):
                        read_pair['discard_reson'] = 'dropped_umi_template'
                    
                    # Check if the UMI has many low quality bases
                    if not read_pair['discard_reson'] and (self.umi_end - self.umi_start) >= self.umi_quality_bases and \
                    len([b for b in read_pair['quality_fw'][self.umi_start:self.umi_end] if (ord(b) - self.phred) < self.min_qual]) > self.umi_quality_bases:
                        read_pair['discard_reson'] = 'dropped_umi'
                
                    # If reverse read has a high AT content discard...
                    if not read_pair['discard_reson'] and self.do_AT_filter and \
                    ((read_pair['sequence_rv'].count("A") + read_pair['sequence_rv'].count("T")) / len(read_pair['sequence_rv'])) * 100 >= self.filter_AT_content:
                        read_pair['discard_reson'] = 'dropped_AT'
                
                    # If reverse read has a high GC content discard...
                    if not read_pair['discard_reson'] and self.do_GC_filter and \
                    ((read_pair['sequence_rv'].count("G") + read_pair['sequence_rv'].count("C")) / len(read_pair['sequence_rv'])) * 100 >= self.filter_GC_content:
                        read_pair['discard_reson'] = 'dropped_GC'
                
                    if not read_pair['discard_reson']:
        
                        if self.do_adaptorA and len(read_pair['sequence_rv']) > self.min_length: read_pair['sequence_rv'], read_pair['quality_rv'] = removeAdaptor(read_pair['sequence_rv'], read_pair['quality_rv'], self.adaptorA, self.adaptor_missmatches) 
                        if self.do_adaptorT and len(read_pair['sequence_rv']) > self.min_length: read_pair['sequence_rv'], read_pair['quality_rv'] = removeAdaptor(read_pair['sequence_rv'], read_pair['quality_rv'], self.adaptorT, self.adaptor_missmatches) 
                        if self.do_adaptorG and len(read_pair['sequence_rv']) > self.min_length: read_pair['sequence_rv'], read_pair['quality_rv'] = removeAdaptor(read_pair['sequence_rv'], read_pair['quality_rv'], self.adaptorG, self.adaptor_missmatches) 
                        if self.do_adaptorC and len(read_pair['sequence_rv']) > self.min_length: read_pair['sequence_rv'], read_pair['quality_rv'] = removeAdaptor(read_pair['sequence_rv'], read_pair['quality_rv'], self.adaptorC, self.adaptor_missmatches)
                        if self.do_adaptorN and len(read_pair['sequence_rv']) > self.min_length: read_pair['sequence_rv'], read_pair['quality_rv'] = removeAdaptor(read_pair['sequence_rv'], read_pair['quality_rv'], self.adaptorN, self.adaptor_missmatches)
                            
                        # Check if the read is smaller than the minimum after removing artifacts   
                        if len(read_pair['sequence_rv']) < self.min_length:
                            read_pair['discard_reson'] = 'dropped_adaptor'
                        else:              
                            # Trim reverse read (will return None if length of trimmed sequence is less than min_length)
                            read_pair['sequence_rv'], read_pair['quality_rv'] = trim_quality(read_pair['sequence_rv'], read_pair['quality_rv'], self.min_qual, self.min_length, self.phred)
                            if not read_pair['sequence_rv'] or not read_pair['quality_rv']:
                                read_pair['discard_reson'] = 'dropped_adaptor'
                    
                    #self.output_read_queue.put( read_pair )
                    out_chunk.append( read_pair )
        
                    count += 1
                    if self.verbose and count % self.stat_line_interwall == 0:
                        self.print_stat_line(start_time, last_time, count, last_count, identity)
                        last_time = time.time()
                        last_count = count
                self.output_read_queue.put( out_chunk )
            #if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: parallel worker process pid='+str(os.getpid())+' - queue empty.\n')

        if self.verbose:
            self.print_stat_line(start_time, last_time, count, last_count, identity)
            last_time = time.time()
            last_count = count

        if self.verbose: sys.stderr.write('InputReadsFilter::INFO:: parallel worker process pid='+str(os.getpid())+' - completed.\n')

    def print_stat_line(self, start_time, last_time, count, last_count, identity, header=False):
        """
        Function that prints a current status row to stderr
        """
        
        import sys
        import time
        
        if header: sys.stderr.write('PROCESS\tREADSPROCESSED\tTIME (s)\tAV_SPEED (reads/s)\tCU_SPEED (reads/s)\n')
        
        sys.stderr.write(
                str(identity)+'\t'+\
                str(count)+'\t'+\
                str(round(time.time()-start_time,3))+'\t'+\
                str(round(count/(time.time()-start_time),2))+'\t'+\
                str(round((count-last_count)/(time.time()-last_time),2))+'\n'
            )