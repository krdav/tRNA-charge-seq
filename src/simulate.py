import copy, os, shutil, bz2, warnings, contextlib, json, random
from subprocess import Popen, PIPE, STDOUT
from mpire import WorkerPool
from Bio import SeqIO
import pandas as pd
import numpy as np
from numpy.random import Generator, PCG64, SFC64, Philox
rng_pg = Generator(PCG64())
import matplotlib.pyplot as plt



class TRNA_ReadSim:
    '''
    Simulate tRNA-Seq reads.
    Initialize by building a tRNA reference set,
    complete with modified positions and shared
    modifications based on sequence similarity.
    Using this reference set reads can be sampled
    by specifying sampling parameters or a sample
    sheet enabling multiple samples to be generated.

    Keyword arguments:
    min_mods -- Minimum modified positions per tRNA (default 5)
    max_mods -- Maximum modified positions per tRNA (default 16)
    dstr_mods -- Distribution to draw the number of modifications from (default {'norm': {'mean': 10, 'std': 3}})
    accp_max_dist -- Maximum hamming distance for sharing modifications i.e. no sharing at this distance (default 10)
    accp_min_dist -- Minimum hamming distance for sharing modifications i.e. certain sharing at this distance (default 2)
    accp_dist_fun -- Function for conversion from hamming distance to probability of sharing modifications (default {'power': {'k': 0.3}})
    ref_len -- Length reference sequences to isolate (default 76)
    ref_excl -- String in name to exclude from reference fx if spike-in control (default 'Escherichia_coli')
    mod_sites -- Positions in the tRNA to modify. If None, using pre-speficied list of 23 positions (default None)
    '''
    def __init__(self, tRNA_ref_fnam, mod_types_fnam, \
                 min_mods=5, max_mods=16, \
                 dstr_mods={'norm': {'mean': 10, 'std': 3}}, \
                 accp_max_dist=10, accp_min_dist=2, \
                 accp_dist_fun={'power': {'k': 0.3}}, \
                 ref_len=76, ref_excl='Escherichia_coli', \
                 mod_sites=None):
        self.tRNA_ref_fnam = tRNA_ref_fnam
        self.ref_len = ref_len
        self.ref_excl = ref_excl
        self.min_mods = min_mods
        self.max_mods = max_mods
        self.dstr_mods = dstr_mods
        self.accp_max_dist = accp_max_dist
        self.accp_min_dist = accp_min_dist
        self.accp_dist_fun = accp_dist_fun
        # Defining a list of possibly modified positions
        # based on where real modifications typically occurs
        # +1 index
        if mod_sites is None:
            self.mod_sites = [6, 8, 9, 10, 16, 17, \
                              20, 22, 26, 27, 32, 34, \
                              35, 37, 39, 40, 46, 47, \
                              48, 49, 54, 55, 58]
        else:
            self.mod_sites = mod_sites
        self.mod_sites = np.array(self.mod_sites)
        self.mod_sites -= 1

        # Define a list of 20 modifications with a probabilities
        # of penetrance and RT events:
        with open(mod_types_fnam, 'r') as fh:
            self.mod_types = eval(fh.read())
        # Make sure the RT event probabilities sum to 1:
        for mod_id in self.mod_types.keys():
            event_prob = np.zeros(10)
            for eid, prob in self.mod_types[mod_id]['RT event'].items():
                event_prob[eid] = prob
            event_prob = event_prob / event_prob.sum()
            self.mod_types[mod_id]['event_ids'] = np.arange(0, 10)
            self.mod_types[mod_id]['event_prob'] = event_prob


        # Sample the reference set:
        self.ref_dct = self._gen_ref()
            
    def _gen_ref(self):
        # Extract all cyto sequences with ref_len nt. length:
        ref_dict = dict()
        for record in SeqIO.parse(self.tRNA_ref_fnam, "fasta"):
            seq = str(record.seq)
            name = str(record.id)
            if len(seq) == self.ref_len and 'mito' not in name and self.ref_excl not in name:
                ref_dict[name] = seq

        # Apply modifications to reference sequences
        # using an iterative approach:
        unmod_pile = copy.deepcopy(ref_dict)
        mod_pile = dict()
        mod_types_arr = np.array(list(self.mod_types.keys()))
        while len(mod_pile) < len(ref_dict):
            for name in list(unmod_pile.keys()):
                if name in mod_pile:
                    continue
                seq = unmod_pile[name]
                del unmod_pile[name]
                mod_pos = self._draw_mod_sites(self.mod_sites, self.min_mods, self.max_mods, self.dstr_mods)
                mod_draw = rng_pg.choice(mod_types_arr, len(mod_pos), replace=True)
                # Make modification array:
                mod_arr = np.zeros(len(seq), dtype=int)
                mod_arr[mod_pos] = mod_draw
                # Make penetrance probability array:
                pen_prob_arr = np.zeros(len(seq), dtype=float)
                pen_prob_arr[mod_pos] = np.array([self.mod_types[m]['penetrance'] for m in mod_draw])
                # Move to pile of modified sequences:
                mod_pile[name] = {'seq': np.array(list(seq)), 'mods': mod_arr, 'pen_prob': pen_prob_arr}

                # Find modification acceptors:
                name_arr = np.array(list(unmod_pile.keys()))
                hdist_arr = np.array([hm_dist(seq, unmod_pile[an]) for an in name_arr])
                prob_arr = self.hamming2prob(hdist=hdist_arr, max_dist=self.accp_max_dist, \
                                             min_dist=self.accp_min_dist, dist_fun=self.accp_dist_fun)
                name_draw = rng_pg.binomial(1, prob_arr) > 0
                # Move those selected to pile of modified sequences:
                for acc_name in name_arr[name_draw]:
                    acc_seq = unmod_pile[acc_name]
                    del unmod_pile[acc_name]
                    # Mask positions where donor nt. != acceptor nt.
                    dnr_mask = np.array([dnt != ant for dnt, ant in zip(seq, acc_seq)])
                    mod_arr_accp = copy.deepcopy(mod_arr)
                    mod_arr_accp[dnr_mask] = 0
                    mod_pile[acc_name] = {'seq': np.array(list(acc_seq)), 'mods': mod_arr_accp, 'pen_prob': pen_prob_arr}

        assert(len(unmod_pile) == 0)
        return(mod_pile)

    def _draw_mod_sites(self, pos_arr, Nmin, Nmax, dist_param):
        if 'norm' in dist_param:
            Nmods = int(round(rng_pg.normal(dist_param['norm']['mean'], dist_param['norm']['std'])))
            if Nmods > Nmax:
                Nmods = Nmax
            elif Nmods < Nmin:
                Nmods = Nmin
        elif 'uniform':
            Nmods = int(round(rng_pg.uniform(Nmin, Nmax)))
        return(np.sort(rng_pg.choice(pos_arr, Nmods, replace=False)))

    def hamming2prob(self, hdist=None, max_dist=None, min_dist=None, \
                     dist_fun=None, plot=False):
        '''
        Convert hamming distance to probability of
        sharing modifications.
        
        Keyword arguments:
        plot -- Only return plot of the resulting distance to probability conversion (default False)
        '''
        if plot is False and any(None is el for el in [hdist, max_dist, min_dist, dist_fun]):
            raise Exception('Must provide all arguments if not plotting: hdist, max_dist, min_dist, dist_fun')
        elif plot:
            pa_lst = list()
            pa_dlst = [self.accp_max_dist, self.accp_min_dist, self.accp_dist_fun]
            for pa_d, pa in zip(pa_dlst, [max_dist, min_dist, dist_fun]):
                if pa is None:
                    pa_lst.append(pa_d)
                else:
                    pa_lst.append(pa)
            max_dist, min_dist, dist_fun = pa_lst

        if 'linear' in dist_fun:
            k = 1
        else:
            k = dist_fun['power']['k']
        h2p = lambda x: (1- (x - min_dist)/(max_dist-min_dist))**k
        if plot:
            hplot = np.arange(0, max_dist+1)
            p = np.ones(max_dist+1)
            p[(min_dist+1):] = h2p(hplot[(min_dist+1):])
            plt.scatter(hplot, p)
            plt.plot(hplot, p)
            plt.xlabel('Hamming distance')
            plt.ylabel('Probability')
            return

        if type(hdist) is int or type(hdist) is float:
            hdist = np.array([hdist])

        hdist_min = hdist < min_dist
        hdist_max = hdist > max_dist
        hdist[hdist_max] = max_dist
        p_arr = h2p(hdist)
        p_arr[hdist_min] = 1.0
        return(p_arr)

    def shuffle_reference(self):
        ''' Re-sample reference database.'''
        self.ref_dct = self._gen_ref()

    def sim_reads(self, btch_size=1e6, bsln_mis=1e-3, \
                     bsln_gap=1e-4, bsln_stop=1e-3, mod_pen_scl=1, \
                     mod_rdthr_scl=0.1, gap_add_lbd=0.4, gap_max=4, \
                     AA_charge=50, strip_gaps=True):
        '''
        Sample tRNA-Seq reads.
        
        Keyword arguments:
        btch_size -- Number of reads to simulate (default 1e6)
        bsln_mis -- Baseline (i.e. not modification induced) mismatch rate (default 1e-3)
        bsln_gap -- Baseline (i.e. not modification induced) gap rate (default 1e-4)
        bsln_stop -- Baseline (i.e. not modification induced) RT stop rate (default 1e-3)
        mod_pen_scl -- Modification penetrance scaling from 0 to 1, with 0 making all modifications non-penetrant(default 1)
        mod_rdthr_scl -- RT readthrough scaling from 0 to 1, with 1 cancelling all RT stop events (default 0.1)
        gap_add_lbd -- Poisson distribution lambda to draw additional gaps downstream a modification induced gap position (default 0.4)
        gap_max -- Max number of gaps induced by one modification (default 4)
        AA_charge -- Aminoacylation charge as defined by the CCA to CC ratio. CCA or CC-ends sampled randomly (default 50)
        strip_gaps -- Strip gap characters '-' from output reads (default True)
        '''

        # Simulate reads using sequence simulation parameters
        # and the two input dictionaries:
        # self.ref_dct => containing arrays of 'seq' (ATGC) and 'mods' (0 if none, 1-20 mod codes)
        # self.mod_types => for each mod code contains penetrance and RT event probabilities
        btch_size = int(btch_size)
        ref_names = list(self.ref_dct.keys())
        ref_idx_arr = np.arange(0, len(ref_names))
        mat_dim = (btch_size, self.ref_len)
        seq_dtype = self.ref_dct[ref_names[0]]['seq'].dtype
        # Storing the reads (i.e. sequences) and modifications
        # in a matrix with one column per read: 
        read_mat = np.empty(mat_dim, dtype=seq_dtype)
        read_anno = list()
        mod_mat = np.zeros(mat_dim, dtype=int)
        # Penetrance matrix:
        pen_prob_mat = np.zeros(mat_dim, dtype=float)

        # Sample reference sequences, insert them into batch matrix:
        ref_idx_draw = rng_pg.choice(ref_idx_arr, btch_size, replace=True)
        for ci, ri in enumerate(ref_idx_draw):
            read_anno.append(ref_names[ri])
            read_mat[ci, :] = self.ref_dct[ref_names[ri]]['seq']
            mod_mat[ci, :] = self.ref_dct[ref_names[ri]]['mods']
            # Scale the penetrance to user input:
            pen_prob_mat[ci, :] = self.ref_dct[ref_names[ri]]['pen_prob'] * mod_pen_scl


        # Draw baseline RT events:
        event_arr = np.array([0, 1, 8, 9]) # 0=no event, 1=random nt., 8=gap, 9=stop
        prob_arr = np.array([0, bsln_mis, bsln_gap, bsln_stop])
        prob_arr[0] = 1 - sum(prob_arr)
        bsln_evnt_mat = rng_pg.choice(event_arr, p=prob_arr, size=mat_dim, replace=True)

        # Sample modified positions based on penetrance:
        mod_mask = rng_pg.binomial(1, pen_prob_mat) > 0

        # Silence non-penetrating modifications:
        mod_mat[~mod_mask] = 0
        # Silence baseline events coinciding with a modification:
        bsln_evnt_mat[mod_mask] = 0

        # Draw RT events on modified positions:
        event_ids = np.arange(0, 10)
        mod_evnt_mat = np.zeros(mat_dim, dtype=int)
        for mt in self.mod_types.keys():
            md_mask = (mod_mat == mt)
            # Scale the readthrough:
            evnt_prob = copy.deepcopy(self.mod_types[mt]['event_prob'])
            evnt_prob[9] *= (1 - mod_rdthr_scl)
            evnt_prob /= evnt_prob.sum()
            mod_evnt_mat[md_mask] = rng_pg.choice(event_ids, p=evnt_prob, size=md_mask.sum())

        # Collect all RT events:
        evnt_mat = bsln_evnt_mat + mod_evnt_mat

        # Sample RT events and insert into each read: 
        rnd_nt = np.array(list('ATGC'))
        pur_nt = np.array(list('AG'))
        pyr_nt = np.array(list('TC'))
        non_nbs = {
            'A': np.array(list('TGC')),
            'T': np.array(list('AGC')),
            'G': np.array(list('ATC')),
            'C': np.array(list('ATG'))
        }
        for rte in event_ids:
            if rte == 0:    # silent
                continue
            elif rte in [5, 6, 7]: # reserved but not used
                continue
            ri, ci = np.where(evnt_mat == rte)
            if rte == 1:    # random nt.
                read_mat[ri, ci] = rng_pg.choice(rnd_nt, size=len(ri))
            elif rte == 2:  # random purine
                read_mat[ri, ci] = rng_pg.choice(pur_nt, size=len(ri))
            elif rte == 3:  # random pyrimidine
                read_mat[ri, ci] = rng_pg.choice(pyr_nt, size=len(ri))
            elif rte == 4:  # non-matching nt.
                for nbs in list('ATGC'):
                    nbs_mask = (read_mat[ri, ci] == nbs)
                    read_mat[ri[nbs_mask], ci[nbs_mask]] = rng_pg.choice(non_nbs[nbs], size=sum(nbs_mask))

            elif rte == 8:  # non-matching nt.
                # Sample number of gaps to expand with:
                gap_nghbr = rng_pg.poisson(gap_add_lbd, size=len(ri))
                gap_nghbr[gap_nghbr > gap_max] = gap_max
                # Expand the gap positions from 3p towards the 5p:
                ri_exp = np.zeros(sum(gap_nghbr) + len(ri), dtype=int)
                ci_exp = np.zeros(sum(gap_nghbr) + len(ri), dtype=int)
                exp_i = 0
                for rii, cii, Ngap in zip(ri, ci, gap_nghbr):
                    for gp in range(Ngap+1):
                        ri_exp[exp_i] = rii
                        ci_exp[exp_i] = cii - gp
                        exp_i += 1
                ci_exp[ci_exp < 0] = 0 # gaps cannot go beyond the 5p (obviously)
                read_mat[ri_exp, ci_exp] = '-'

            elif rte == 9:  # stop
                # Find the maximum column index
                # for each stop i.e. the stop
                # closest to the 3p:
                max_ci_dict = dict()
                for rii, cii in zip(ri, ci):
                    try:
                        max_ci_dict[rii] = max(max_ci_dict[rii], cii)
                    except:
                         max_ci_dict[rii] = cii

                # Enforce the stop by removing bases in the read:
                for rii, cii in max_ci_dict.items():
                    read_mat[rii, :(cii+1)] = ''

        # Trim 3p base off all uncharged tRNA reads:
        unch_prob = (100 - AA_charge) / 100
        unch_mask = rng_pg.binomial(1, unch_prob, size=mat_dim[0]) > 0
        read_mat[unch_mask, -1] = ''
        
        # Return simulated reads with their reference annotations:
        if strip_gaps:
            reads = [''.join(read_mat[si]).replace('-', '') for si in range(mat_dim[0])]
        else:
            # Still strip 5p gaps but preserve internal gaps:
            reads = [''.join(read_mat[si]).lstrip('-') for si in range(mat_dim[0])]
        return(read_anno, reads)

    def adjust_ref_rdthr(self, rdthr_max=40, rdthr_min=30, N_btch=5, \
                         btch_size=1e4, max_try=50, strip_gaps=False,
                         kwargs={}):
        '''
        Re-sample reference set until
        achieving the requested readthrough.
        
        Keyword arguments:
        rdthr_max --  Maximum percentage readthrough in the N_btch samples (default 43)
        rdthr_min --  Minimum percentage readthrough in the N_btch samples (default 37)
        N_btch --  Number of read batches to sample to generate min/max readthrough (default 5)
        btch_size -- Number of reads to simulate to generate min/max readthrough (default 1e4)
        max_try --  Maximum number of iteration before giving up (default 50)
        strip_gaps -- Strip gap characters '-' from output reads (default False)
        '''

        attempt = 0
        while attempt <= max_try:
            self.shuffle_reference()
            cov_lst = list()
            for _ in range(N_btch):
                ann, reads = self.sim_reads(btch_size=btch_size, strip_gaps=strip_gaps, 
                                            AA_charge=100, **kwargs)
                cov = 100 * sum(len(sr) == self.ref_len for sr in reads)/len(reads)
                cov_lst.append(cov)
            if max(cov_lst) <= rdthr_max and min(cov_lst) >= rdthr_min:
                print('Stopped with reference having {}/{} max/min coverage accross {} read batches.'.format(max(cov_lst), min(cov_lst), N_btch))
                break
            attempt += 1
        if attempt > max_try:
            print('Failed to adjust reference to requested readthrough')

    def _handle_row(self, index, row):
        '''Handle the read simulation for each sample row.'''

        # Find the read simulation parameters set
        # in the simulation sample sheet:
        params = {kwargs: row[kwargs] for kwargs in self.sim_seq_kwargs if kwargs in row}
        # Generate the simulated reads:
        ann_lst, reads = self.sim_reads(**params)
        # Generate UMIs:
        if not self.UMI_ns:
            Nseqs = len(reads)
            umi_mat = np.empty((Nseqs, 10), dtype=np.dtype('U1'))
            chN = np.array(list('ATGC'), dtype=np.dtype('U1'))
            chPyr = np.array(list('TC'), dtype=np.dtype('U1'))
            umi_mat[:, 0:9] = rng_pg.choice(chN, size=(Nseqs, 9))
            umi_mat[:, 9:10] = rng_pg.choice(chPyr, size=(Nseqs, 1))
            umi_lst = [''.join(ur) for ur in umi_mat]
        else:
            Nseqs = len(reads)
            chN = np.array(list('ATGC'), dtype=np.dtype('U1'))
            umi_mat1 = rng_pg.choice(chN, size=(Nseqs, 1))
            umi_mat2 = rng_pg.choice(chN, size=(Nseqs, 2))
            umi_mat3 = rng_pg.choice(chN, size=(Nseqs, 3))
            umi_lst = [''.join(random.choice([u1, u2, u3])) for u1, u2, u3 in zip(umi_mat1, umi_mat2, umi_mat3)]

        # Write simulated reads in fastq format:
        anno_json = dict()
        fastq_fnam = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, row['sample_name_unique'])
        with bz2.open(fastq_fnam, "wt") as fastq_fh:
            idx = 0
            for umi, ann, seq in zip(umi_lst, ann_lst, reads):
                read_id = '{}_{}'.format(row['sample_name_unique'], idx)
                anno_json[read_id] = ann
                # Add barcode and UMI sequence to title:
                title = '{} {}:{}'.format(read_id, row['barcode_seq'], umi)
                # Add Phred score basically
                # meanining no sequencing error:
                qual = 'J' * len(seq)
                # Write the simulated sequence in fastq format:
                fastq_fh.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))
                idx += 1

        # Write annotations for each simulated read:
        anno_fnam = '{}/{}.json.bz2'.format(self.anno_dir_abs, row['sample_name_unique'])
        with bz2.open(anno_fnam, "wt") as fh_out:
            json.dump(anno_json, fh_out)
        
        # Return row with UMI stats:
        row['N_after_downsample'] = Nseqs
        row['N_UMI_observed'] = len(set(umi_lst))
        UMI_bins = 4**9*2
        row['N_UMI_expected'] = UMI_bins*(1-((UMI_bins-1) / UMI_bins)**Nseqs)
        row['percent_UMI_obs-vs-exp'] = 100 * row['N_UMI_observed']/row['N_UMI_expected']
        return(row)

    def sim_from_sheet(self, sim_sheet_fnam, NBdir, n_jobs=4, \
                       data_dir='sim_data', overwrite_dir=False,
                       UMI_ns=False):
        '''
        Simulate reads with difference parameters
        using a simulation sample sheet.
        
        Keyword arguments:
        n_jobs -- Number of multiprocessing jobs to use (default 4)
        data_dir -- Name of folder to put simulated reads in. If not existing it will be created (default 'sim_data')
        UMI_ns -- Lazy way of implementing a different UMI specification
        '''

        # Keyword arguments for the read simulation parameters
        # to search for in the simulation sample sheet:
        self.UMI_ns = UMI_ns
        self.sim_seq_kwargs = ['btch_size', 'bsln_mis', 'bsln_gap', \
                               'bsln_stop', 'mod_pen_scl', 'mod_rdthr_scl', \
                               'gap_add_lbd', 'gap_max', 'AA_charge', 'strip_gaps']

        # Make folder structure:
        self.dir_dict = {'NBdir': NBdir, 'data_dir': data_dir, 'UMI_dir': 'UMI_trimmed', 'anno_dir': 'read_anno'}
        # Make data dir:
        data_dir_abs = '{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'])
        if not os.path.exists(data_dir_abs):
            os.mkdir(data_dir_abs)

        # Make UMI dir:
        self.UMI_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['UMI_dir'])
        try:
            os.mkdir(self.UMI_dir_abs)
        except:
            if overwrite_dir:
                shutil.rmtree(self.UMI_dir_abs)
                os.mkdir(self.UMI_dir_abs)
            else:
                raise Exception('UMI_trimmed already exists under the folder {} and overwrite is set to False'.format(self.dir_dict['data_dir']))
        
        # Make simulation annotation dir:
        self.anno_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['anno_dir'])
        try:
            os.mkdir(self.anno_dir_abs)
        except:
            if overwrite_dir:
                shutil.rmtree(self.anno_dir_abs)
                os.mkdir(self.anno_dir_abs)
            else:
                raise Exception('read_anno already exists under the folder {} and overwrite is set to False'.format(self.dir_dict['data_dir']))

        # Generate simulated reads for each sample
        # in the simulation sample sheet:
        sim_df = pd.read_excel(sim_sheet_fnam)
        data = list(sim_df.iterrows())
        '''
        Persistent strange behaviour when using parallel 
        processing with random sampling. Maybe related
        to this:
        https://github.com/numpy/numpy/issues/9650

        def initfn():
            np.random.seed()
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self._handle_row, data, worker_init=initfn)
        '''
        results = [self._handle_row(index, row) for index, row in sim_df.iterrows()]

        # Make and write new sample sheet:
        sample_df_fnam_abs = '{}/new_sample_list.xlsx'.format(self.dir_dict['NBdir'])
        sample_df = pd.DataFrame(results)
        sample_df = sample_df.drop(columns=self.sim_seq_kwargs, errors='ignore')
        sample_df.to_excel(sample_df_fnam_abs, index=False)
        sp_set = set(sample_df['species'])
        assert(len(sp_set) == 1)
        species = sp_set.pop()

        # Dump new reference:
        tRNA_db_sim = self.write_sim_tRNA_db(self.dir_dict, species)

        return(self.dir_dict, sample_df, tRNA_db_sim)

    def write_sim_tRNA_db(self, dir_dict, species, out_dir='tRNA_database_sim'):
        '''
        Write the tRNA database used for simulation with folder structure,
        sequences and BLAST index such that it can be used
        directly for an alignment run.
        '''

        # Make out_dir folder:
        out_dir_abs = '{}/{}/{}'.format(dir_dict['NBdir'], dir_dict['data_dir'], out_dir)
        try:
            os.mkdir(out_dir_abs)
        except:
            shutil.rmtree(out_dir_abs)
            os.mkdir(out_dir_abs)

        tRNA_db_sim = dict()
        # Make out_dir species folder:
        out_dir_sp_abs = '{}/{}'.format(out_dir_abs, species)
        try:
            os.mkdir(out_dir_sp_abs)
        except:
            shutil.rmtree(out_dir_sp_abs)
            os.mkdir(out_dir_sp_abs)
        
        # Write the masked sequences as fasta:
        out_fnam = '{}-tRNAs.fa'.format(species)
        out_fnam_abs = '{}/{}'.format(out_dir_sp_abs, out_fnam)
        tRNA_db_sim[species] = out_fnam_abs
        with open(out_fnam_abs, 'w') as fh:
            for name in self.ref_dct.keys():
                seq = ''.join(self.ref_dct[name]['seq'])
                fh.write('>{}\n{}\n'.format(name, seq))

        # Make BLAST DB on fasta:
        self._makeblastdb(out_dir_sp_abs, out_fnam, self.dir_dict['NBdir'])

        return(tRNA_db_sim)

    def _makeblastdb(self, file_dir, fnam, NBdir):
        '''
        Build BLAST database (required for SWIPE).
        Notice, it is built with sequence type = protein
        This is to enable a custom scoring matrix for SWIPE.
        '''
        try:
            os.chdir(file_dir)
            cmd_str = 'makeblastdb -dbtype prot -in {} -blastdb_version 4'.format(fnam)
            cmd = cmd_str.split()
            log_fn = 'makeblastdb_logfile.txt'
            with Popen(cmd, stdout=PIPE, stderr=STDOUT) as p, open(log_fn, 'a') as file:
                file.write('Starting subprocess with command:')
                file.write(str(cmd))
                file.write('\n')
                for line in p.stdout:
                    file.write(line.decode('utf-8'))
                file.write('\n****** DONE ******\n\n\n')
            
            os.chdir(NBdir)
        except Exception as err:
            os.chdir(NBdir)
            raise err







try:
    import jellyfish
    def hm_dist(s1, s2):
        return(jellyfish.hamming_distance(s1, s2))
except:
    def hm_dist(s1, s2):
        return(sum(c1!=c2 for c1, c2 in zip(s1, s2)))