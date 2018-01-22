import sys
import os
import logging
import multiprocessing
from collections import defaultdict
import numpy as np
import pandas as pd
import pysam
from stpipeline.common.clustering import *
from stpipeline.common.unique_events_parser import uniqueEventsParser
from stpipeline.common.dataset import computeUniqueUMIs
from stpipeline.common.gff_reader import gff_lines
from stpipeline.common.sam_utils import get_annotations

class DatasetCreator():

    def __init__(self,
                  input_file,
                  qa_stats,
                  gff_filename,
                  umi_cluster_algorithm="hierarchical",
                  umi_allowed_mismatches=1,
                  umi_counting_offset=250,
                  output_folder=None,
                  output_template=None,
                  verbose=True):
        """
        The functions parses the reads in BAM format
        that have been annotated and demultiplexed (containing spatial barcode).
        It then groups them by gene-barcode to count reads accounting for duplicates
        using the UMIs (clustering them suing the strand and start position).
        It outputs the records in a matrix of counts in TSV format and BED format and it also
        writes out some statistics.
        :param input_file: the file with the annotated-demultiplexed records in BAM format
        :param qa_stats: the Stats object to add some stats (THIS IS PASSED BY REFERENCE)
        :param umi_cluster_algorithm: the clustering algorithm to cluster UMIs
        :param umi_allowed_mismatches: the number of miss matches allowed to remove
                                      duplicates by UMIs
        :param umi_counting_offset: the number of bases allowed as offset (start position) when counting UMIs
        :param output_folder: path to place the output files
        :param output_template: the name of the dataset
        :param verbose: True if we can to collect the stats in the logger
        :type input_file: str
        :type umi_cluster_algorithm: str
        :type umi_allowed_mismatches: boolean
        :type umi_counting_offset: integer
        :type output_folder: str
        :type output_template: str
        :type verbose: bool
        :raises: RuntimeError,ValueError,OSError,CalledProcessError
        """

        self.input_file = input_file
        self.qa_stats = qa_stats
        self.gff_filename = gff_filename
        self.umi_cluster_algorithm = umi_cluster_algorithm
        self.umi_allowed_mismatches = umi_allowed_mismatches
        self.umi_counting_offset = umi_counting_offset
        self.output_folder = output_folder
        self.output_template = output_template
        self.verbose = verbose

        self.logger = logging.getLogger("STPipeline")

        if not os.path.isfile(self.input_file):
            error = "Error creating dataset, input file not present {}\n".format(self.input_file)
            self.logger.error(error)
            raise RuntimeError(error)

        if self.output_template:
            self.filenameDataFrame = "{}_stdata.tsv".format(self.output_template)
            self.filenameReadsBED = "{}_reads.bed".format(self.output_template)
        else:
            self.filenameDataFrame = "stdata.tsv"
            self.filenameReadsBED = "reads.bed"

        # Obtain the clustering function
        if self.umi_cluster_algorithm == "naive":
            self.group_umi_func = countUMINaive
        elif self.umi_cluster_algorithm == "hierarchical":
            self.group_umi_func = countUMIHierarchical
        elif self.umi_cluster_algorithm == "Adjacent":
            self.group_umi_func = dedup_adj
        elif self.umi_cluster_algorithm == "AdjacentBi":
            self.group_umi_func = dedup_dir_adj
        elif self.umi_cluster_algorithm == "Affinity":
            self.group_umi_func = affinity_umi_removal
        else:
            error = "Error creating dataset.\n" \
            "Incorrect clustering algorithm {}".format(self.umi_cluster_algorithm)
            self.logger.error(error)
            raise RuntimeError(error)

    def run(self, ):
        print 'starting run'

        print 'create workers'
        threads = multiprocessing.cpu_count()
        self.workers = [ worker_process(
                                self.group_umi_func,
                                self.umi_allowed_mismatches,
                                self.umi_counting_offset,
                                self.output_folder,
                                verbose=self.verbose
                            ) for i in range(threads) ]

        print 'create gene controller'
        self.gene_controller = gene_controller(self.gff_filename, [self.input_file], self.qa_stats, verbose=self.verbose)
        self.gene_controller.get_genes()
        self.gene_controller.connect_to_workers(self.workers)
        worker_pids=[worker.process.pid for worker in self.workers]
        print worker_pids

        counts_table = self.gene_controller.collect_results()

        # Write data frame to file
        counts_table.to_csv(os.path.join(self.output_folder, self.filenameDataFrame), sep="\t", na_rep=0)

        with open(os.path.join(self.output_folder, self.filenameReadsBED), "w") as reads_handler:
            for pid in worker_pids:
                infile = open(os.path.join(self.output_folder,'reads.worker_{}.bed'.format( pid )))
                reads_handler.write( infile.read() )
                infile.close()
                os.remove( infile.name )

        # connect them by queues/pipes
        # start them
        # monior their progress
        # kill them
        # clean up

class gene_controller():

    def __init__(self, gff_filename, input_file_names, qa_stats, verbose=False):
        # create queues and other things
        self.gff_filename = gff_filename
        self.input_file_names = input_file_names
        self.genes = dict()
        self.workers = []
        self.results_queue = multiprocessing.Queue()
        self.logger = logging.getLogger("STPipeline")
        self.qa_stats = qa_stats
        self.verbose=verbose

    def get_genes(self,):

        # Read gtf and get genes
        self.get_gene_coordinates(self.gff_filename)

        # fetch read count per gene from bam files
            # also get ambiguos gene annotations from bam files
            # genes = { 'gene_id':{'name':gene_name, 'id':gene_id, 'sequence':chrom, 'start':start, 'end':end, 'read_count':read_count, 'files':{filename:read_count}} ...}
        self.get_read_counts()

    def get_gene_coordinates(self, gff_filename):
        """
        function that reads the coordinates of genes
        and save them as values of a dictionary with the gene ID as key
        dict: [GENE_id] => gene_coordinate
        """
        print 'get gene coordinates from gtf'
        cdef dict genes
        cdef dict line
        cdef str gene_id
        cdef str seqname
        cdef int start
        cdef int end

        genes = dict()

        for line in gff_lines(gff_filename):

            seqname = line['seqname']
            start = int(line['start'])
            end = int(line['end'])
            gene_id = line['gene_id']

            try:
                if gene_id[0] == '"' and gene_id[-1] == '"': gene_id=gene_id[1:-1]
            except KeyError:
                raise ValueError(
                    'The gene_id attribute is missing in the annotation file ({0})\n'.format(self.gff_filename)
                    )
            try:
                if int(start) < genes[ gene_id ]['start']:genes[ gene_id ][ 'start' ]= int(end)
                if int(end)   > genes[ gene_id ]['end']:  genes[ gene_id ][ 'end' ]  = int(end)
            except KeyError:
                genes[ gene_id ] = {
                    'id':gene_id,
                    'sequence':seqname,
                    'start':start,
                    'end':end,
                    'read_count':0,
                    'files':{}
                    }

        self.genes = genes

    def get_read_counts(self,):

        print 'fetch annotations from bam files'
        keyword_arguments = { 'return_queue':multiprocessing.Queue(), 'return_filename':True }
        sub_processes = [
                multiprocessing.Process(target=get_annotations, args=[bam_file_name], kwargs=keyword_arguments)
                for bam_file_name in self.input_file_names
            ]

        for process in sub_processes: process.start()

        for bam_file_name in self.input_file_names: pysam.index( bam_file_name, '{0}.bai'.format(bam_file_name) )

        def empty_queue():
            while not keyword_arguments['return_queue'].empty():
                temp_dict,filename = keyword_arguments['return_queue'].get()
                for annotation,count in temp_dict.iteritems():
                    if annotation == '__no_feature': continue # THIS SHOULD BE FIXED SOMEHOW!! OR IGNORED BUT NOTED
                    try:
                        assert filename not in self.genes[annotation]['files']
                        self.genes[annotation]['read_count'] += count
                        self.genes[annotation]['files'][filename] = count
                    except KeyError:
                        if annotation[0:len('__ambiguous[')] == '__ambiguous[':#'GENE1+GENE2]'
                            try: # to get the right most right most coordinate ;)
                                ambiguous_gene_ids = annotation[len('__ambiguous['):-1].split('+')
                                sequence_names = [ self.genes[gene_id]['sequence'] for gene_id in ambiguous_gene_ids ]
                                assert len(set(sequence_names)) == 1
                                self.genes[ annotation ] = {
                                        'id':annotation,
                                        'sequence':sequence_names[0],
                                        'start':min([ self.genes[gene_id]['start'] for gene_id in ambiguous_gene_ids ]),
                                        'end':max([ self.genes[gene_id]['end'] for gene_id in ambiguous_gene_ids ]),
                                        'read_count':count,
                                        'files':{filename:count}
                                        }
                            except KeyError:
                                raise ValueError('ERROR:: gene with id {0} is not found in gtf file\n'.format(annotation))
                        else:
                            raise ValueError('ERROR:: gene with id {0} is not found in gtf file\n'.format(annotation))

                        # IMPORTANT!!!! we do not capture any "__no_feature" annotations is this a problem?

        while True in [process.is_alive() for process in sub_processes]: empty_queue()

        for process in sub_processes: process.join()

        empty_queue()

    def connect_to_workers(self, workers):
        # sort genes by readcount and distribute to workers
            # put genes in queue for workers to process

        for worker in workers:
            worker.genes=[]
            worker.return_queue = self.results_queue

        i = 0
        key=lambda x: x['read_count']
        for gene in sorted(self.genes.values(),key=key,reverse=True):
            workers[i].genes.append(gene)
            i +=1
            if i == len(workers): i=0

        for worker in workers:
            worker.process.start()

    def empty_queue(self, ):
        while not self.results_queue.empty():
            gene_id, transcript_counts_by_spot, discarded_reads_worker = self.results_queue.get()
            print 'got gene {} in results collection'.format( gene_id )
            # Add spot and dict [gene] -> count to containers
            self.list_indexes.append(gene_id)
            self.list_row_values.append(transcript_counts_by_spot)
            self.total_record += sum([transcript_count for spot_coordinate,transcript_count in transcript_counts_by_spot.iteritems() ])
            self.discarded_reads += discarded_reads_worker

    def collect_results(self, ):

        # Containers needed to create the data frame
        self.list_row_values = list()
        self.list_indexes = list()
        # Some counters
        self.total_record = 0
        self.discarded_reads = 0

        # get gene spot combos from worker
        i=1
        while True in [worker.process.is_alive() for worker in self.workers]:
            print 'empty queue call #',i;
            i +=1
            self.empty_queue()
        for worker in self.workers: worker.process.join()
        print 'empty queue call #',i,'(final)';
        self.empty_queue()

        # Print spot gene combos to tsv file (mem or file)
        # Create the data frame
        counts_table = pd.DataFrame(self.list_row_values, index=self.list_indexes)
        counts_table.fillna(0, inplace=True)
        counts_table=counts_table.T # Transpose the dictionary to still get the spots as rows and genes as columns in the final tsv
        #print counts_table

        # Compute some statistics
        total_record = self.total_record
        discarded_reads = self.discarded_reads
        total_barcodes = len(counts_table.index)
        total_transcripts = np.sum(counts_table.values, dtype=np.int32)
        number_genes = len(counts_table.columns)
        max_count = counts_table.values.max()
        min_count = counts_table.values.min()
        aggregated_spot_counts = counts_table.sum(axis=1).values
        aggregated_gene_counts = (counts_table != 0).sum(axis=1).values
        max_genes_feature = aggregated_gene_counts.max()
        min_genes_feature = aggregated_gene_counts.min()
        max_reads_feature = aggregated_spot_counts.max()
        min_reads_feature = aggregated_spot_counts.min()
        average_reads_feature = np.mean(aggregated_spot_counts)
        average_genes_feature = np.mean(aggregated_gene_counts)
        std_reads_feature = np.std(aggregated_spot_counts)
        std_genes_feature = np.std(aggregated_gene_counts)

        # Print some statistics
        if self.verbose:
            self.logger.info("Number of unique molecules present: {}".format(total_transcripts))
            self.logger.info("Number of unique events (gene-feature) present: {}".format(total_record))
            self.logger.info("Number of unique genes present: {}".format(number_genes))
            self.logger.info("Max number of genes over all features: {}".format(max_genes_feature))
            self.logger.info("Min number of genes over all features: {}".format(min_genes_feature))
            self.logger.info("Max number of unique molecules over all features: {}".format(max_reads_feature))
            self.logger.info("Min number of unique molecules over all features: {}".format(min_reads_feature))
            self.logger.info("Average number genes per feature: {}".format(average_genes_feature))
            self.logger.info("Average number unique molecules per feature: {}".format(average_reads_feature))
            self.logger.info("Std number genes per feature: {}".format(std_genes_feature))
            self.logger.info("Std number unique molecules per feature: {}".format(std_reads_feature))
            self.logger.info("Max number of unique molecules over all unique events: {}".format(max_count))
            self.logger.info("Min number of unique molecules over all unique events: {}".format(min_count))
            self.logger.info("Number of discarded reads (possible duplicates): {}".format(discarded_reads))

        # Update the QA object
        self.qa_stats.reads_after_duplicates_removal = int(total_transcripts)
        self.qa_stats.unique_events = total_record
        self.qa_stats.barcodes_found = total_barcodes
        self.qa_stats.genes_found = number_genes
        self.qa_stats.duplicates_found = discarded_reads
        self.qa_stats.max_genes_feature = max_genes_feature
        self.qa_stats.min_genes_feature = min_genes_feature
        self.qa_stats.max_reads_feature = max_reads_feature
        self.qa_stats.min_reads_feature = min_reads_feature
        self.qa_stats.max_reads_unique_event = max_count
        self.qa_stats.min_reads_unique_event = min_count
        self.qa_stats.average_gene_feature = average_genes_feature
        self.qa_stats.average_reads_feature = average_reads_feature

        return counts_table

class worker_process():

    def __init__(self, group_umi_func, umi_allowed_mismatches, umi_counting_offset, output_folder, verbose=False):

        self.genes = None
        self.process = multiprocessing.Process(target=self.worker_loop)
        self.output_folder = output_folder
        self.group_umi_func = group_umi_func
        self.umi_allowed_mismatches = umi_allowed_mismatches
        self.umi_counting_offset = umi_counting_offset
        self.logger = logging.getLogger("STPipeline")
        self.verbose=verbose

    def get_gene_data(self, gene):

        #print gene['id']
        if gene['read_count'] == 0: return {} #return gene['id'], []

        # define variables cython style for more efficient looping and in memory storage
        cdef object rec
        cdef str clear_name
        cdef int mapping_quality
        cdef int start
        cdef int end
        cdef str chrom
        cdef str strand
        cdef tuple transcript
        cdef tuple spot_coordinates
        cdef int x
        cdef int y
        cdef str gene_id
        cdef dict spots = {}

        reads_total = 0
        for filename, count in gene['files'].iteritems():

            # get the genes from the bamfiles
            # pysam.index( filename, '{0}.bai'.format(filename) )
            sam_file = pysam.AlignmentFile(filename, "rb")
            reads_in_bam = 0

            # Parse the coordinate sorted bamfile record by record i.e. by genome
            # coordinate from first chromosome to last
            for rec in sam_file.fetch( gene['sequence'], gene['start']-100, gene['end']+100 ):

                # get the info about the record from the bam file
                clear_name = rec.query_name
                mapping_quality = rec.mapping_quality

                # Account for soft-clipped bases when retrieving the start/end coordinates
                start = int(rec.reference_start - rec.query_alignment_start)
                end = int(rec.reference_end + (rec.query_length - rec.query_alignment_end))
                chrom = sam_file.getrname(rec.reference_id)
                strand = "+"
                if rec.is_reverse:
                    # We swap start and end if the transcript mapped to the reverse strand
                    strand = "-"
                    start, end = end, start

                # Get TAGGD tags from the bam file
                x,y,gene_id,umi = (-1,-1,'None','None')
                for (k, v) in rec.tags:
                    if k == "B1":
                        x = int(v) ## The X coordinate
                    elif k == "B2":
                        y = int(v) ## The Y coordinate
                    elif k == "XF":
                        gene_id = str(v) ## The gene name
                    elif k == "B3":
                        umi = str(v) ## The UMI
                    else:
                        continue
                # Check that all tags are present
                if 'None' in [gene_id,umi] or -1 in [x,y]:
                    logger.warning("Warning parsing annotated reads.\n" \
                                   "Missing attributes for record {}\n".format(clear_name))
                    continue

                if gene['id'] != gene_id: continue

                # Create a new transcript and add it to the in memory gene_buffer dictionary
                transcript = (chrom, start, end, clear_name, mapping_quality, strand, umi)
                spot_coordinates = (x,y)

                try:
                    spots[spot_coordinates].append(transcript)
                except KeyError:
                    spots[spot_coordinates] = [transcript]

                reads_in_bam += 1

            # Close the bam file and yield the last gene(s)
            sam_file.close()

            assert reads_in_bam == count, 'ASSERTION ERROR:: {}, reads_in_bam {} != {} count \n'.format( gene_id, reads_in_bam, count )
            reads_total += reads_in_bam

        assert reads_total == gene['read_count']

        #print 'returning',gene['id']
        #print 'spots=',spots
        return spots

    @property
    def pid(self,):
        return self.process.pid

    def worker_loop(self,):
        print 'kicking off worker {}'.format( os.getpid() )
        # do the umi clustering etc for one gene
        # put results back to controller
        # Parse unique events to generate the unique counts and the BED file
        #unique_events_parser = uniqueEventsParser(self.input_file, self.gff_filename)
        filenameReadsBED = 'reads.worker_{}.bed'.format(os.getpid())
        total_record = 0
        with open(os.path.join(self.output_folder, filenameReadsBED), "w") as reads_handler:
            # this is the generator returning a dictionary with spots for each gene
            #for gene, spots in unique_events_parser.all_unique_events():
            for gene in self.genes:
                spots = self.get_gene_data( gene )
                #print 'gene=',gene, 'spot=',spots
                if not spots: continue

                transcript_counts_by_spot = {}
                discarded_reads = 0
                for spot_coordinates, reads in spots.iteritems():
                    (x,y) = spot_coordinates
                    # Re-compute the read count accounting for duplicates using the UMIs
                    # Transcripts is the list of transcripts (chrom, start, end, clear_name, mapping_quality, strand, UMI)
                    # First:
                    # Get the original number of transcripts (reads)
                    read_count = len(reads)
                    # Compute unique transcripts (based on UMI, strand and start position +- threshold)
                    unique_transcripts = computeUniqueUMIs(reads, self.umi_counting_offset,
                                                           self.umi_allowed_mismatches, self.group_umi_func)
                    # The new transcript count
                    transcript_count = len(unique_transcripts)
                    assert transcript_count > 0 and transcript_count <= read_count
                    # Update the discarded reads count
                    discarded_reads += (read_count - transcript_count)
                    # Update read counts in the container (replace the list
                    # of transcripts for a number so it can be exported as a data frame)
                    transcript_counts_by_spot["{0}x{1}".format(x, y)] = transcript_count
                    # Write every unique transcript to the BED output (adding spot coordinate and gene name)
                    for read in unique_transcripts:
                        reads_handler.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(read[0],
                                                                                                   read[1],
                                                                                                   read[2],
                                                                                                   read[3],
                                                                                                   read[4],
                                                                                                   read[5],
                                                                                                   gene['id'],
                                                                                                   x,y))
                    # keep a counter of the number of unique events (spot - gene) processed
                    total_record += 1

                self.return_queue.put( (gene['id'], transcript_counts_by_spot, discarded_reads) )
                #print 'worker {}: gene {} put in queue'.format( os.getpid(), gene['id'] )

        if total_record == 0:
            error = "Error creating dataset, input file did not contain any transcript\n"
            self.logger.error(error)
            raise RuntimeError(error)