import sys, os, shutil, bz2, json
from subprocess import Popen, PIPE, STDOUT
import xml.etree.ElementTree as ET
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import numpy as np
from mpire import WorkerPool
from json_stream import streamable_dict
import json_stream



class SWIPE_align:
    '''
    This class is used to align reads to a reference database.
    In the case of tRNAs, a fasta file with known tRNA sequences
    for the appropriate species.
    This method uses the SWIPE aligner:
    https://github.com/torognes/swipe
    SWIPE performs Smith-Waterman local sequence alignment
    on each read against all sequences in the database i.e.
    the results do not rely on any alignemnt heuristics
    which other aligners like BLAST, Bowtie, GSNAP etc do.
    Because no heuristics are used, the results are guaranteed
    to contain the best alignment, as defined by the alignment score.
    '''
    def __init__(self, dir_dict, tRNA_database, sample_df, score_mat, gap_penalty=6, extension_penalty=1, min_score_align=20, common_seqs=None):
        # Swipe command template: 
        self.swipe_cmd_tmp = 'swipe\t--query\tINPUT_FILE\t--db\tDATABASE_FILE\t--out\tOUTPUT_FILE\t--symtype\t1\t--outfmt\t7\t--num_alignments\t3\t--num_descriptions\t3\t--evalue\t0.000000001\t--num_threads\t12\t--strand\t1\t--matrix\tSCORE_MATRIX\t-G\tGAP_PENALTY\t-E\tEXTENSION_PENALTY'
        self.swipe_cmd_tmp = self.swipe_cmd_tmp.replace('SCORE_MATRIX', score_mat)
        self.swipe_cmd_tmp = self.swipe_cmd_tmp.replace('GAP_PENALTY', str(gap_penalty))
        self.swipe_cmd_tmp = self.swipe_cmd_tmp.replace('EXTENSION_PENALTY', str(extension_penalty))
        self.SWIPE_overwrite = True
        self.dry_run = False
        # Input:
        self.tRNA_database, self.sample_df, self.min_score_align = tRNA_database, sample_df, min_score_align
        self.dir_dict = dir_dict
        self.common_seqs_fnam = common_seqs
        self.common_seqs_dict = dict() # map common sequences to index

        # Check files exists before starting:
        self.UMI_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['UMI_dir'])
        for _, row in self.sample_df.iterrows():
            trimmed_fn = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, row['sample_name_unique'])
            assert(os.path.exists(trimmed_fn))

        # If using common sequences, check the input format and existence,
        # then read into dictionary:
        if self.common_seqs_fnam is not None:
            # We can only allow one species if using common sequences.
            # Multiple species would require running the alignment on common sequences
            # several times, defeating the purpose, but also making the code much
            # more complicated.
            sp_set = set(self.sample_df['species'].values)
            if len(sp_set) > 1:
                raise Exception('Only one species allowed in sample sheet when using common sequences.')
            self.common_seqs_sp = list(sp_set)[0]

            print('Using common sequences to prevent duplicated alignment.')
            assert(os.path.exists(self.common_seqs_fnam))
            assert(self.common_seqs_fnam[-4:] == '.bz2')
            with bz2.open(self.common_seqs_fnam, "rt") as input_fh:
                for ridx, record in enumerate(SeqIO.parse(input_fh, "fasta")):
                    seq = str(record.seq)
                    assert(ridx == int(record.id))
                    assert(not seq in self.common_seqs_dict)
                    self.common_seqs_dict[seq] = ridx
            self.Ncommon = len(self.common_seqs_dict)

    def make_dir(self, overwrite=True):
        # Create folder for files:
        self.align_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['align_dir'])
        try:
            os.mkdir(self.align_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.align_dir_abs)
                os.mkdir(self.align_dir_abs)
            else:
                print('Using existing folder because overwrite set to false: {}'.format(self.align_dir_abs))
        return(self.align_dir_abs)

    def run_parallel(self, n_jobs=4, overwrite=True, verbose=True):
        self.SWIPE_overwrite = overwrite
        self.verbose = verbose
        if self.verbose:
            print('Running Swipe on:', end='')
        os.chdir(self.align_dir_abs)
        try:
            # Run SWIPE in parallel:
            data = list(self.sample_df.iterrows())
            if not self.common_seqs_fnam is None:
                data.append((0, 'common-seqs'))
                self.__prep_common()
            with WorkerPool(n_jobs=n_jobs) as pool:
                swipe_return = pool.map(self.__start_SWIPE, data)
            if not self.common_seqs_fnam is None:
                os.remove(self.common_seqs_fnam[:-4])
            # Collect results in parallel:
            if self.verbose:
                print('\nCollecting alignment statistics, from sample:', end='')
            with WorkerPool(n_jobs=n_jobs) as pool:
                results = pool.map(self.__collect_stats, data)
            self.__write_stats(results)
            os.chdir(self.dir_dict['NBdir'])
            return(self.sample_df)
        except Exception as err:
            os.chdir(self.dir_dict['NBdir'])
            raise err
    
    '''
    Currently, I do not want to support serial run, instead use n_jobs=1.
    Keeping this legacy code because it can be usefull for error handling.

    def run_serial(self, dry_run=False, overwrite=True, verbose=True):
        self.dry_run = dry_run
        self.SWIPE_overwrite = overwrite
        self.verbose = verbose
        if self.verbose:
            print('Running Swipe on:', end='')
        os.chdir(self.align_dir_abs)
        try:
            swipe_return = [self.__start_SWIPE(index, row) for index, row in self.sample_df.iterrows()]
            # Collect results:
            if not dry_run:
                if self.verbose:
                    print('\nCollecting alignment statistics, from sample:', end='')
                results = [self.__collect_stats(index, row) for index, row in self.sample_df.iterrows()]
                self.__write_stats(results)
            os.chdir(self.dir_dict['NBdir'])
            return(self.sample_df)
        except Exception as err:
            os.chdir(self.dir_dict['NBdir'])
            raise err
    '''

    def __prep_common(self):
        # Convert reads to uncompressed fasta as required by Swipe:
        with bz2.open(self.common_seqs_fnam, 'rb') as fh_in:
            with open(self.common_seqs_fnam[:-4], 'wb') as fh_out:
                shutil.copyfileobj(fh_in, fh_out)

    def __start_SWIPE(self, index, row):
        # Generate Swipe command line:
        if type(row) == str:  # special case for common sequences
            sample_name_unique = 'common-seqs'
            sp_tRNA_database = self.tRNA_database[self.common_seqs_sp]
            swipe_cmd, swipe_outfile = self.__make_SWIPE_cmd(sp_tRNA_database, self.common_seqs_fnam[:-4], sample_name_unique)
        else:
            sp_tRNA_database = self.tRNA_database[row['species']]
            sample_name_unique = row['sample_name_unique']
            trimmed_fn = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, sample_name_unique)
            trimmed_fasta_fn = trimmed_fn[:-10] + '.fasta'
            swipe_cmd, swipe_outfile = self.__make_SWIPE_cmd(sp_tRNA_database, trimmed_fasta_fn, sample_name_unique)
        if self.dry_run:
            print('Swipe cmd: {}'.format(' '.join(swipe_cmd)))
            return(1)

        # Skip, if results file has already been made and no overwrite:
        SWres_fnam = '{}_SWalign.json.bz2'.format(sample_name_unique)
        SWnohits_fnam = '{}_SWalign-nohits.fasta.bz2'.format(sample_name_unique)
        if not self.SWIPE_overwrite and os.path.isfile(SWres_fnam) and os.path.isfile(SWnohits_fnam):
            return(1)

        # Prepare sequences for SWIPE:
        if type(row) == str:
            pass # common sequences have already been prepared
        elif not self.common_seqs_fnam is None:
            # Count the number of times a common sequence is observed:
            common_obs = np.zeros(self.Ncommon)
            with bz2.open(trimmed_fn, 'rt') as fh_in:
                with open(trimmed_fasta_fn, 'wt') as fh_out:
                    for title, seq, qual in FastqGeneralIterator(fh_in):
                        if seq in self.common_seqs_dict: # Count commont sequence
                            seq_idx = self.common_seqs_dict[seq]
                            common_obs[seq_idx] += 1
                        else: # Write sequence not found in common
                            fh_out.write('>{}\n{}\n'.format(title, seq))
            # Write the common sequence observations to a file:
            common_obs_fn = '{}_common-seq-obs.json'.format(sample_name_unique)
            with open(common_obs_fn, 'w') as fh_out:
                json.dump(common_obs.tolist(), fh_out)
        else: # do not use common sequences
            # Convert reads to fasta as required by Swipe:
            with bz2.open(trimmed_fn, 'rt') as fh_bz:
                SeqIO.convert(fh_bz, "fastq", trimmed_fasta_fn, 'fasta')

        # Submit the process:
        if self.verbose:
            print('  {}'.format(sample_name_unique), end='')
        log_fn = '{}_logfile.txt'.format(sample_name_unique)
        with Popen(swipe_cmd, stdout=PIPE, stderr=STDOUT) as p, open(log_fn, 'a') as file:
            file.write('Starting subprocess with command:')
            file.write(str(swipe_cmd))
            file.write('\n')
            for line in p.stdout:
                file.write(line.decode('utf-8'))
            file.write('\n****** DONE ******\n\n\n')

        # Reformat output so it can be read as XML:
        swipe_outfile_xml = self.__prep_SWIPE_XML(swipe_outfile)
        # Read XML file as a streamable generator:
        json_stream = self.__parse_SWIPE_XML(swipe_outfile_xml, sp_tRNA_database)
        # Dump query_hits as JSON:
        with bz2.open(SWres_fnam, 'wt', encoding="utf-8") as fh:
            json.dump(json_stream, fh)

        # Remove tmp files:
        if type(row) != str:
            os.remove(trimmed_fasta_fn)
        os.remove(swipe_outfile)
        os.remove(swipe_outfile_xml)
        return(1)
    
    def __make_SWIPE_cmd(self, sp_tRNA_database, trimmed_fn, sample_name_unique):
        swipe_cmd = self.swipe_cmd_tmp
        swipe_cmd = swipe_cmd.replace('DATABASE_FILE', sp_tRNA_database)
        swipe_cmd = swipe_cmd.replace('INPUT_FILE', trimmed_fn)
        swipe_outfile = '{}_SWalign'.format(sample_name_unique)
        swipe_cmd = swipe_cmd.replace('OUTPUT_FILE', swipe_outfile)
        swipe_cmd = swipe_cmd.split('\t')
        return(swipe_cmd, swipe_outfile)
    
    def __prep_SWIPE_XML(self, swipe_outfile):
        # Add "data" as root for the xml file:
        swipe_outfile_xml = swipe_outfile + '.xml'
        xml_first_line = '<data>\n'
        xml_last_line = '</data>\n'
        with open(swipe_outfile, 'r') as from_file:
            try:
                os.remove(swipe_outfile_xml)
            except:
                pass
            with open(swipe_outfile_xml, 'a') as to_file:
                from_file.readline()
                to_file.write(xml_first_line)
                shutil.copyfileobj(from_file, to_file)
                to_file.write(xml_last_line)
        return(swipe_outfile_xml)

    # Use decorator to define this as a streamable dict
    # to avoid loading all into memory when writing to JSON.
    # See: https://pypi.org/project/json-stream/
    @streamable_dict
    def __parse_SWIPE_XML(self, swipe_outfile_xml, sp_tRNA_database):
        # Read the database IDs and use them to verify alignment results:
        db_id_set = set()
        for record in SeqIO.parse(sp_tRNA_database, "fasta"):
            db_id_set.add(record.id)
        
        # Parse XML:
        hit_dict = {tag: [] for tag in ['score', 'query', 'name', 'qpos', 'dpos', 'qseq', 'aseq', 'dseq']}
        pickup = True # When True, pick up hit data and store in tmp dict ("hit_dict")
        flush = False # When True, flush tmp dict into "query_hits"
        high_score = -999
        hit_dict_prev = None # For debugging
        SWxml = ET.iterparse(swipe_outfile_xml)

        for event, elem in SWxml:
            # When "result" tag is encountered it marks the end of the hits for a query.
            # Flush the data picked up:
            if elem.tag == 'result':
                # This would not save any memory:
                # elem.clear() # clear for saving memory
                flush = True
            # Pick up all tags defined in "hit_dict":
            elif pickup and elem.tag in hit_dict:
                hit_dict[elem.tag].append(elem.text)
                # If all highest alignment score(s) have been picked up,
                # stop picking up more data:
                if elem.tag == 'score':
                    if int(elem.text) >= high_score:
                        high_score = int(elem.text)
                    else:
                        pickup = False

            # Flush out hit results into "query_hits".
            # Only if results are stored and alignment score is above minimum:
            if flush and len(hit_dict['score']) > 0 and high_score >= self.min_score_align:
                query_hits = dict()
                # Convert alignment score to integers:
                hit_dict['score'] = [int(s) for s in hit_dict['score']]
                # Find all the highest scoring hits, extract indices for selection:
                high_score_idx = indices(hit_dict['score'], high_score)
                # Remove all hits with alignment score lower than
                # the maximun score:
                for tag in hit_dict:
                    hit_dict[tag] = [hit_dict[tag][hidx] for hidx in high_score_idx]
                # Convert qpos/dpos (query/database alignment position) string to integer tuple:
                hit_dict['qpos'] = [tuple(map(int, qp.split(','))) for qp in hit_dict['qpos']]
                hit_dict['dpos'] = [tuple(map(int, dp.split(','))) for dp in hit_dict['dpos']]
                # Assert that only one query sequence has been picked up:
                ls_query = list(set(hit_dict['query']))
                assert(len(ls_query) == 1)
                # Start to populate the dict entry for the query sequence:
                query = ls_query[0]
                query_hits['score'] = high_score
                # The "name" tag is the database result.
                # First extract the right hand side of the string,
                # corresponding to the fasta header,
                # then sort (if multiple hits) and merge with @:
                hit_dict['name'] = [n.split(' ')[-1] for n in hit_dict['name']]
                for n in hit_dict['name']: # quick assertion that name is in database
                     assert(n in db_id_set)
                # Extract sorting index for other data to be sorted:
                name_idx = sorted(range(len(hit_dict['name'])), key=lambda k: hit_dict['name'][k])
                name = '@'.join([hit_dict['name'][didx] for didx in name_idx])
                query_hits['name'] = name
                # Add qpos/dpos:
                query_hits['qpos'] = [hit_dict['qpos'][didx] for didx in name_idx]
                query_hits['dpos'] = [hit_dict['dpos'][didx] for didx in name_idx]
                # Add alignment strings, but only for the first hit:
                query_hits['qseq'] = hit_dict['qseq'][name_idx[0]]
                query_hits['aseq'] = hit_dict['aseq'][name_idx[0]]
                query_hits['dseq'] = hit_dict['dseq'][name_idx[0]]
                query_hits['aligned'] = True

                yield query, query_hits
            elif flush:
                ls_query = list(set(hit_dict['query']))
                if len(ls_query) > 0:
                    query_hits = dict()
                    query = ls_query[0]
                    query_hits['aligned'] = False
                    yield query, query_hits

            # After flushing, reset the variables for new data pickup:
            if flush:
                hit_dict_prev = hit_dict.copy() # For debugging
                hit_dict = {tag: [] for tag in ['score', 'query', 'name', 'qpos', 'dpos', 'qseq', 'aseq', 'dseq']}
                flush = False
                pickup = True
                high_score = -999

    def __collect_stats(self,  index, row):
        # Collect stats about the alignment #
        if type(row) == str:
            sample_name_unique = 'common-seqs'
        else:
            sample_name_unique = row['sample_name_unique']
        # If common sequences were used get their observations:
        if type(row) != str and not self.common_seqs_fnam is None:
            common_obs_fn = '{}_common-seq-obs.json'.format(sample_name_unique)
            with open(common_obs_fn, 'r') as fh_in:
                common_obs = json.load(fh_in)

        if self.verbose:
            print('  {}'.format(sample_name_unique), end='')

        # Collect information:
        query_nohits = set()
        N_mapped = 0
        N_mult_mapped = 0
        # Read query_hits from JSON:
        SWres_fnam = '{}_SWalign.json.bz2'.format(sample_name_unique)
        with bz2.open(SWres_fnam, 'rt', encoding="utf-8") as SWres_fh:
            # Parse JSON data as a stream,
            # i.e. as a transient dict-like object
            SWres = json_stream.load(SWres_fh)
            for readID, align_dict in SWres.persistent().items():
                if not align_dict['aligned']:
                    query_nohits.add(readID)
                else:
                    N_mapped += 1
                    if '@' in align_dict['name']:
                        N_mult_mapped += 1

        # Collect information from common sequences:
        if type(row) != str and not self.common_seqs_fnam is None:
            # Read query_hits from JSON:
            SWres_fnam = '{}_SWalign.json.bz2'.format('common-seqs')
            with bz2.open(SWres_fnam, 'rt', encoding="utf-8") as SWres_fh:
                # Parse JSON data as a stream,
                # i.e. as a transient dict-like object
                SWres = json_stream.load(SWres_fh)
                for readID, align_dict in SWres.persistent().items():
                    if align_dict['aligned']:
                        readID_int = int(readID)
                        N_mapped += common_obs[readID_int]
                        if '@' in align_dict['name']:
                            N_mult_mapped += common_obs[readID_int]

        # Dump unaligned sequences:
        SWnohits_fnam = '{}_SWalign-nohits.fasta.bz2'.format(sample_name_unique)
        if type(row) == str:
            with bz2.open(SWnohits_fnam, 'wt', encoding="utf-8") as fh_out:
                with bz2.open(self.common_seqs_fnam, 'rt') as fh_in:
                    for record in SeqIO.parse(fh_in, "fasta"):
                        seq_id = record.id
                        if seq_id in query_nohits:
                            fh_out.write(">{}\n{}\n".format(record.id, str(record.seq)))
            return(False)
        else:
            trimmed_fn = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, sample_name_unique)
            with bz2.open(SWnohits_fnam, 'wt', encoding="utf-8") as fh_out:
                with bz2.open(trimmed_fn, 'rt') as fh_in:
                    for title, seq, qual in FastqGeneralIterator(fh_in):
                        seq_id = title.split()[0]
                        if seq_id in query_nohits:
                            fh_out.write(">{}\n{}\n".format(title, seq))

            # Calculate stats:
            if row['N_after_trim'] == 0:
                map_p = 0
            else:
                map_p = N_mapped / row['N_after_trim'] * 100
            # Multiple mappings have fasta IDs merged with "@":
            if N_mapped == 0:
                P_ma = 0
            else:
                P_ma = N_mult_mapped / N_mapped * 100
            P_sa = 100 - P_ma

            return([sample_name_unique, N_mapped, P_sa, P_ma, map_p])

    def __write_stats(self, results):
        # Remove the results entry from common seqeunces:
        results = [res for res in results if not res is False]
        stats_df = pd.DataFrame(results, columns=['sample_name_unique', 'N_mapped', 'percent_single_annotation', 'percent_multiple_annotation', 'Mapping_percent'])
        # Merge stats with sample info dataframe:
        self.sample_df = self.sample_df.drop(columns=['N_mapped', 'percent_single_annotation', 'percent_multiple_annotation', 'Mapping_percent'], errors='ignore')
        self.sample_df = self.sample_df.merge(stats_df, on=['sample_name_unique'])
        # Write stats:
        self.sample_df.to_excel('sample_stats.xlsx')


def indices(lst, element):
    result = []
    offset = -1
    while True:
        try:
            offset = lst.index(element, offset+1)
        except ValueError:
            return result
        result.append(offset)

