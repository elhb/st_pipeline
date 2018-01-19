import sys
import os
import logging
import multiprocessing
from collections import defaultdict
import numpy as np
import pandas as pd
from stpipeline.common.clustering import *
from stpipeline.common.unique_events_parser import uniqueEventsParser
from stpipeline.common.dataset import computeUniqueUMIs
from stpipeline.common.gff_reader import gff_lines

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
            filenameDataFrame = "{}_stdata.tsv".format(self.output_template)
            filenameReadsBED = "{}_reads.bed".format(self.output_template)
        else:
            filenameDataFrame = "stdata.tsv"
            filenameReadsBED = "reads.bed"
             
        # Some counters
        total_record = 0
        discarded_reads = 0
        
        # Obtain the clustering function
        if self.umi_cluster_algorithm == "naive":
            group_umi_func = countUMINaive
        elif self.umi_cluster_algorithm == "hierarchical":
            group_umi_func = countUMIHierarchical
        elif self.umi_cluster_algorithm == "Adjacent":
            group_umi_func = dedup_adj
        elif self.umi_cluster_algorithm == "AdjacentBi":
            group_umi_func = dedup_dir_adj
        elif self.umi_cluster_algorithm == "Affinity":
            group_umi_func = affinity_umi_removal
        else:
            error = "Error creating dataset.\n" \
            "Incorrect clustering algorithm {}".format(self.umi_cluster_algorithm)
            self.logger.error(error)
            raise RuntimeError(error)
     
        # Containers needed to create the data frame
        list_row_values = list()
        list_indexes = list()   
    
        # Parse unique events to generate the unique counts and the BED file    
        unique_events_parser = uniqueEventsParser(self.input_file, self.gff_filename)
        with open(os.path.join(self.output_folder, filenameReadsBED), "w") as reads_handler: ######################################## DO WE EVEN NEEED THE BEDFILE!?!?! ###############
            # this is the generator returning a dictionary with spots for each gene
            for gene, spots in unique_events_parser.all_unique_events(): 
                transcript_counts_by_spot = {}
                for spot_coordinates, reads in spots.iteritems():
                    (x,y) = spot_coordinates
                    # Re-compute the read count accounting for duplicates using the UMIs
                    # Transcripts is the list of transcripts (chrom, start, end, clear_name, mapping_quality, strand, UMI)
                    # First:
                    # Get the original number of transcripts (reads)
                    read_count = len(reads)
                    # Compute unique transcripts (based on UMI, strand and start position +- threshold)
                    unique_transcripts = computeUniqueUMIs(reads, self.umi_counting_offset, 
                                                           self.umi_allowed_mismatches, group_umi_func)
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
                                                                                                   gene,
                                                                                                   x,y)) 
                    # keep a counter of the number of unique events (spot - gene) processed
                    total_record += 1
                    
                # Add spot and dict [gene] -> count to containers
                list_indexes.append(gene)
                list_row_values.append(transcript_counts_by_spot)
                
        if total_record == 0:
            error = "Error creating dataset, input file did not contain any transcript\n"
            self.logger.error(error)
            raise RuntimeError(error)
        
        # Create the data frame
        counts_table = pd.DataFrame(list_row_values, index=list_indexes)
        counts_table.fillna(0, inplace=True)
        counts_table=counts_table.T # Transpose the dictionary to still get the spots as rows and genes as columns in the final tsv
        
        # Compute some statistics
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
         
        # Write data frame to file
        counts_table.to_csv(os.path.join(self.output_folder, filenameDataFrame), sep="\t", na_rep=0)       


class main_controller():

    def __init__(self, threads, gff_filename, input_file_names):
        self.threads = threads
        self.gff_filename = gff_filename
        self.input_file_names = input_file_names

    def run(self, ):
        
        # create workers
        self.workers = [worker_process() for i in range(threads)]
        
        # create gene controller
        self.gene_controller = gene_controller(self.gff_filename, self.input_file_names)
        self.gene_controller.get_genes()
        self.gene_controller.get_read_counts()
        self.gene_controller.connect_to_workers()
        
        counts_table = self.gene_controller.collect_results()
        
        # connect them by queues/pipes
        # start them
        # monior their progress
        # kill them
        # clean up
        pass

class gene_controller():

    def __init__(self, gff_filename, input_file_names):
        # create queues and other things
        self.gff_filename = gff_filename
        self.input_file_names = input_file_names
        self.genes = dict()
        self.workers = []
        self.results_queue = multiprocessing.Queue()

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
        # from stpipeline.common.sam_utils import get_annotations
        get_annotations = lambda ('gene_id',100,'pelle.bam') # mock function

        keyword_arguments = { 'return_queue':multiprocessing.Queue(), 'return_filename'=True }
        sub_processes = [
                multiprocessing.Process(target=get_annotations, args=[bam_file_name], kwargs=keyword_arguments)
                for bam_file_name in files_to_merge
            ]

        for process in sub_processes: process.start()

        def empty_queue():
            while not keyword_arguments['return_queue'].empty():
                temp_dict = keyword_arguments['return_queue'].get()
                for annotation,count,filename in temp_dict.iteritems():
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
                                self.gene[ annotation ] = {
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

        key=lambda x: x['read_count']
        for gene in sorted(self.genes.values(),key=key,reverse=True)]:
            workers[i].genes.append(gene)
            i +=1
            if i == len(workers): i=0
        
        for worker in workers:
            worker.start()

    def collect_results(self, ):
        # get gene spot combos from worker
        # Print spot gene combos to tsv file (mem or file)
        pass

class worker_process():

    def __init__(self,):
        self.genes = None
        self._process = multiprocessing.Process(target=self.worker_loop)

    def get_genes(self,):
        # get the genes from the bamfiles
        for gene in genes:
            yield gene

    def worker_loop(self,):
        # do the umi clustering etc for one gene
        # put results back to controller
        for gene in self.get_genes():
            print gene


