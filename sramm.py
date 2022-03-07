# !/usr/bin/env python

"""
SRAMM: Short Read Alignment Mapping Metrics
"""

__author__ = "Alvin Chon"
__email__ = "achon@iastate.edu"
__version__ = "0.1.0"

import argparse
import time
import numpy as np
import pysam
from math import floor
from multiprocessing import Pool

# global
output_prefix = ""


def parse_vg_file(vg_file, ofh, end, num_reads, num_cpu):
    """
    Must be the VG output file from the -v option and look like this:
    read_name   seq_name    pos mapq    alignment_score
    ERR174324.2079  chr7    63592360        0       111
    ERR174324.2079  chr5    39196216        0       60
    ERR174324.2079          0       0       0
    :param vg_file: vg debug output
    :param ofh: output file handle
    :param end: str, single of paired
    :param num_reads: number of reads to process as a batch
    :param num_cpu: number of CPUs to use in the pool
    :return: data, obs_max
    """
    in_file = open(vg_file, 'r')
    ctime = time.time()
    pool = Pool(processes=num_cpu)
    chunk_size = int(num_reads / num_cpu / 10)
    prev_read_id = ""
    obs_max = 0
    counter_reads = 0
    data = {}
    mapq_data = {}
    if end == "single":
        for line in in_file:
            line_parts = line.strip().split('\t')
            if len(line_parts) == 5:
                cur_read_id, seq_name, pos, mapq, score_aln = line_parts
            else:
                cur_read_id, seq_name, pos, mapq, score_aln = line_parts[0], "", 0, 0, 0
            score_aln = int(score_aln)
            if score_aln > obs_max:
                obs_max = score_aln
            if cur_read_id == prev_read_id:
                data[cur_read_id] = data[cur_read_id] + [score_aln]
                mapq_data[cur_read_id] = mapq_data[cur_read_id] + [mapq]
            else:
                counter_reads += 1
                if counter_reads >= num_reads:
                    # print('---Starting batch process: ', np.around(time.time() - ctime, 2))
                    batch_single_calc(pool, chunk_size, ofh, data, mapq_data, obs_max, end)
                    data.clear()
                    mapq_data.clear()
                    counter_reads = 0
                data[cur_read_id] = [score_aln]
                mapq_data[cur_read_id] = [mapq]
                prev_read_id = cur_read_id
        print('---Starting final read processing: ', np.around(time.time() - ctime, 2))
        batch_single_calc(pool, chunk_size, ofh, data, mapq_data, obs_max, end)
        pool.close()
        pool.join()
    elif end == "paired":
        paired_index = -1
        paired_scores = -1
        paired_mapq = -1
        for line in in_file:
            line_parts = line.strip().split('\t')
            if len(line_parts) == 5:
                cur_read_id, seq_name, pos, mapq, score_aln = line_parts
            else:
                cur_read_id, seq_name, pos, mapq, score_aln = line_parts[0], "", 0, 0, 0
            score_aln = int(score_aln)
            if score_aln > obs_max:
                obs_max = score_aln
            if cur_read_id == prev_read_id:
                if paired_index == 0:
                    paired_scores = score_aln
                    paired_mapq = mapq
                    paired_index = 1
                elif paired_index == 1:
                    data[cur_read_id] = data[cur_read_id] + [(paired_scores, score_aln)]
                    mapq_data[cur_read_id] = mapq_data[cur_read_id] + [paired_mapq]    # VG has both left and right alns in a pair share the same MAPQ
                    paired_scores = -1
                    paired_mapq = -1
                    paired_index = 0
                    # calc_scores((cur_read_id, data[cur_read_id], obs_max, end_type))
            else:
                counter_reads += 1
                if counter_reads >= num_reads:
                    # print('---Starting batch process: ', np.around(time.time() - ctime, 2))
                    batch_paired_vg_calc(pool, chunk_size, ofh, data, mapq_data, obs_max, end)
                    data.clear()
                    mapq_data.clear()
                    counter_reads = 0
                paired_scores = score_aln
                paired_mapq = mapq
                paired_index = 1
                prev_read_id = cur_read_id
                data[cur_read_id] = []
                mapq_data[cur_read_id] = []
        print('---Starting final read processing: ', np.around(time.time() - ctime, 2))
        batch_paired_vg_calc(pool, chunk_size, ofh, data, mapq_data, obs_max, end)
        pool.close()
        pool.join()
    return obs_max


def parse_alignment_file(input_file, ofh, end, num_reads, num_cpu):
    """
    Parses sam/bam/cram file.  File must be samtools sorted by QNAME using the -n option or with collate.
    :param input_file: bam file
    :param ofh: output file handle
    :param end: str, single or paired
    :param num_reads: number of reads to batch process
    :param num_cpu: number of cpu for MP
    :return:
    """
    ctime = time.time()
    pool = Pool(processes=num_cpu)
    chunk_size = int(num_reads / num_cpu / 10)
    batch = {}
    aln_data = {}  # paired only
    prev_read_id = ""
    obs_max = 0
    counter_reads = 0
    alignment_file = pysam.AlignmentFile(input_file, 'rb')
    if end == 'single':
        data = {}
        mapq_data = {}
        for read in alignment_file:
            score_aln = int(read.get_tag('ASl'))
            cur_read_id = read.query_name
            if score_aln > obs_max:
                obs_max = score_aln
            if cur_read_id == prev_read_id:
                data[cur_read_id] = data[cur_read_id] + [score_aln]
                mapq_data[cur_read_id] = mapq_data[cur_read_id] + [read.mapping_quality]
            else:
                counter_reads += 1
                if counter_reads >= num_reads:
                    # print('---Starting batch process: ', np.around(time.time() - ctime, 2))
                    batch_single_calc(pool, chunk_size, ofh, data, mapq_data, obs_max, end)
                    data.clear()
                    mapq_data.clear()
                    counter_reads = 0
                data[cur_read_id] = [score_aln]
                mapq_data[cur_read_id] = [read.mapping_quality]
                prev_read_id = cur_read_id
        print('---Starting final read processing: ', np.around(time.time() - ctime, 2))
        batch_single_calc(pool, chunk_size, ofh, data, mapq_data, obs_max, end)
        pool.close()
        pool.join()
    elif end == 'paired':
        for read in alignment_file:
            score_aln = int(read.get_tag('AS'))
            cur_read_id = read.query_name
            if score_aln > obs_max:
                obs_max = score_aln
            if cur_read_id == prev_read_id:
                aln_data['{}_{}'.format(read.reference_id, read.reference_start)] = (
                    read.is_paired, read.is_proper_pair, score_aln,
                    read.reference_id, read.reference_start,
                    read.next_reference_id, read.next_reference_start, read.mapping_quality)
            else:
                counter_reads += 1
                if counter_reads >= num_reads:
                    # print('---Starting batch process: ', np.around(time.time() - ctime, 2))
                    batch_paired_calc(pool, chunk_size, ofh, batch, obs_max, end)
                    batch.clear()
                    counter_reads = 0
                else:
                    if prev_read_id != "":  # initial check
                        batch[prev_read_id] = aln_data
                # reset aln_data for new read
                aln_data = {'{}_{}'.format(read.reference_id, read.reference_start): (
                    read.is_paired, read.is_proper_pair,
                    score_aln, read.reference_id,
                    read.reference_start,
                    read.next_reference_id,
                    read.next_reference_start, read.mapping_quality)}  # resets the dict to remove previous read's data
                prev_read_id = cur_read_id
        print('---Starting final read processing: ', np.around(time.time() - ctime, 2))
        batch_paired_calc(pool, chunk_size, ofh, batch, obs_max, end)
        pool.close()
        pool.join()
    return


def batch_single_calc(pool, chunk_size, ofh, data, mapq_data, obs_max, end):
    """
    MP wrapper
    """
    for r in pool.imap_unordered(calc_scores, [(k, data[k], obs_max, end) for k in data.keys()],
                                 chunksize=chunk_size):
        k, map_score, uniq_rep_score, unmap_score = r
        write_stats_line((k, data[k], mapq_data[k], (map_score, uniq_rep_score, unmap_score)), end, ofh)
    return


def batch_paired_vg_calc(pool, chunk_size, ofh, data, mapq_data, obs_max, end):
    """
    MP wrapper
    """
    for r in pool.imap_unordered(process_pairs_vg, [(k, v, mapq_data[k], obs_max, end) for k, v in data.items()],
                                 chunksize=chunk_size):
        write_stats_line(r, end, ofh)
    return


def batch_paired_calc(pool, chunk_size, ofh, batch, obs_max, end):
    """
    MP wrapper
    """
    for r in pool.imap_unordered(process_pairs, [(k, v, obs_max, end) for k, v in batch.items()],
                                 chunksize=chunk_size):
        write_stats_line(r, end, ofh)
    return


def calc_scores(scores_tuple):
    """
    Wrapper for calculating scores
    :param scores_tuple: Tuple containing key, scores, obs_max, end
    :return: read_id and metric scores
    """
    read_id, alignment_scores, obs_max, end = scores_tuple
    map_score = calc_mapping_score(alignment_scores, obs_max, end)
    uniq_rep_score = calc_unique_repeat_score(alignment_scores, obs_max, end)
    unmap_score = calc_unmapped_score(alignment_scores, obs_max, end)
    return read_id, map_score, uniq_rep_score, unmap_score


def process_pairs_vg(alns_tuple):
    """
    Finds the pairs in alignments of one read
    :param alns_tuple: alignments of one read in as a tuple
    :return: read_id, AS pairs, mapq, metric scores
    """
    read_id, alignments, mapqs, obs_max, end = alns_tuple
    return read_id, alignments, mapqs, calc_scores((read_id, alignments, obs_max, end))


def process_pairs(alns_tuple):
    """
    Finds the pairs in alignments of one read
    :param alns_tuple: alignments of one read in as a tuple
    :return: read_id, AS pairs, mapq, metric scores
    """
    read_id, alignments, obs_max, end = alns_tuple
    pairs = []
    mapq = []
    already_paired = []
    for k, aln in alignments.items():
        # print(k, aln)
        paired, proper_pair, score_aln, start_id, start, next_start_id, next_start, mapping_quality = aln
        mate_aln = alignments['{}_{}'.format(next_start_id, next_start)]
        mate_paired, mate_proper_pair, mate_score_aln, mate_start_id, mate_start, mate_next_start_id, mate_next_start, mate_mapping_quality = mate_aln
        if (mate_aln, aln) not in already_paired:
            already_paired.append((aln, mate_aln))
            pairs += [(score_aln, mate_score_aln)]
            mapq.append(mapping_quality)
    return read_id, pairs, mapq, calc_scores((read_id, pairs, obs_max, end))


def calc_mapping_score(alignments, obs_max, end):
    """
    Returns the max alignment score.
    :param alignments: list, alignments for a given read
    :param obs_max: int, obs_max alignment score
    :param end: str, single of paired
    :return: score
    """
    if end == "single":
        return float(max(alignments)) / obs_max
    if end == "paired":
        pair_sum = -1
        for pair in alignments:
            if sum(pair) > pair_sum:
                pair_sum = sum(pair)
        return float(pair_sum) / (2 * obs_max)


def calc_unique_repeat_score(alignments, obs_max, end):
    """
    Returns the unique-repeat score
    :param alignments: list, alignments for a given read
    :param obs_max: int, obs_max alignment score
    :param end: str, single of paired
    :return: score
    """
    if end == "single":
        alignments = sorted(alignments, reverse=True)
        num_alns = len(alignments)
        if num_alns == 0:
            return 0.0
        elif num_alns >= 1:
            base = float(alignments[0]) / (obs_max * num_alns)
            res = 0.0
            for i in range(num_alns):
                res += obs_max - alignments[i]
            res /= obs_max * num_alns
            return base + res
        return -1
    if end == "paired":
        alignments = sorted(alignments, key=lambda s: s[0] + s[1], reverse=True)
        num_alns = len(alignments)
        if num_alns == 0:
            return 0.0
        elif num_alns >= 1:
            base = float(sum(alignments[0])) / (obs_max * 2 * num_alns)
            res = 0.0
            for i in range(num_alns):
                res += obs_max * 2 - sum(alignments[i])
            res /= obs_max * 2 * num_alns
            return base + res
        return -1


def calc_unmapped_score(alignments, obs_max, end):
    """
    Returns the unmapped score
    :param alignments: list, alignments for a given read
    :param obs_max: int, obs_max alignment score
    :param end: str, single of paired
    :return: score
    """
    if end == "single":
        alignments = sorted(alignments, reverse=True)
        num_alns = len(alignments)
        if num_alns == 0:
            return 1.0
        else:
            return 1 - float(sum(alignments)) / (obs_max * num_alns)
    if end == "paired":
        alignments = sorted(alignments, key=lambda s: s[0] + s[1], reverse=True)
        num_alns = len(alignments)
        if num_alns == 0:
            return 1.0
        else:
            return 1 - float(sum([score[0] + score[1] for score in alignments])) / (obs_max * 2 * num_alns)


def write_stats_line(line, end, ofh):
    """
    Writes out the stats file after generation
    :param line: info for one read
    :param end: end type
    :param ofh: output file handle
    :return:
    """
    if end == "single":
        read_id, alns, mapq, stats = line
        try:
            ofh.write("{}\t{}\t{}\t{}\t{}\n".format(read_id, ','.join(str(d) for d in alns),
                                                    len(alns), ','.join(str(m) for m in mapq),
                                                    ','.join(str(np.around(r, 5)) for r in stats)))
        except:
            print('Failed writing {}.  Line information {}'.format(read_id, line))
    if end == "paired":
        read_id, pairs, mapq, stats = line
        stats = stats[1:]  # due to returning read_id for single end
        try:
            ofh.write("{}\t{}\t{}\t{}\t{}\n".format(read_id, ','.join(str(d) for d in pairs),
                                                    len(pairs), ','.join(str(m) for m in mapq),
                                                    ','.join(str(np.around(r, 5)) for r in stats)))
        except:
            print('Failed writing {}.  Line information {}'.format(read_id, line))
    return


def filter_stats_wrapper(num_reads, num_cpu, score_ranges, stat_file, ofh):
    """
    Filters the reads based upon the metric scores and number of alignments
    :param num_reads: number of reads per batch
    :param num_cpu: num of cpus to be used in MP pool
    :param score_range: 4 ranges for MI, UR, UM, and NA
    :param stat_file: file containing the SRAMM metrics generated from stats
    :param ofh: output file handle
    :return:
    """
    pool = Pool(processes=num_cpu)
    chunk_size = int(num_reads / num_cpu)
    ifh = open(stat_file, 'r')
    lines = []
    counter = 0
    for line in ifh:
        lines.append(line)
        counter += 1
        if counter >= num_reads:
            batch_filter(pool, chunk_size, ofh, lines, score_ranges)
            lines = []
            counter = 0
    batch_filter(pool, chunk_size, ofh, lines, score_ranges)
    pool.close()
    pool.join()
    return


def batch_filter(pool, chunk_size, ofh, lines, score_ranges):
    """
    MP wrapper
    """
    for r in pool.imap_unordered(filter_stats, [(line, score_ranges) for line in lines], chunk_size):
        if r is not None:
            ofh.write(r)
    return


def filter_stats(tup):
    """
    Filter reads based upon metrics
    :param tup: contains a line from the stats file and the score ranges to filter by
    :return: line if passing, none otherwise
    """
    line, score_ranges = tup
    read_id, aln_scores, num_alns, mapq_scores, sramm_scores = line.strip().split('\t')
    map_score, uniq_rep_score, unmap_score = [float(v) for v in sramm_scores.split(',')]
    map_range, uniq_rep_range, unmap_range, num_alns_range = score_ranges
    # Tempted to use all() here, but might be more complex looking since ranges are dynamic
    if map_range[0] <= map_score <= map_range[1] and \
            uniq_rep_range[0] <= uniq_rep_score <= uniq_rep_range[1] and \
            unmap_range[0] <= unmap_score <= unmap_range[1] and \
            num_alns_range[0] <= int(num_alns) <= num_alns_range[1]:
        return line
    else:
        return None


def read_stat_file(stat_file, num_reads, num_cpu, mapq_max, num_score_bins):
    """
    Read sramm output stat file for graph generation
    Histogram batching and aggregation to reduce memory size
    :param stat_file: sramm stat file
    :param num_reads: number of reads per batch
    :param num_cpu: num of cpus to be used in MP pool
    :param mapq_max: max value of mapq for graph purposes
    :param num_score_bins: number of bins for the 3 SRAMM metrics.
    :return:
    """
    pool = Pool(processes=num_cpu)
    chunk_size = int(num_reads / num_cpu)
    ifh = open(stat_file, 'r')
    stats = {}
    mapq_hist = np.zeros(mapq_max)  # MAPQ is an 8 bit flag hence 0 to 255 at max...
    mapped_hist = np.zeros(num_score_bins)
    uniq_rep_hist = np.zeros(num_score_bins)
    unmapped_hist = np.zeros(num_score_bins)
    counter = 0
    for line in ifh:
        read_id, aln_scores, num_alns, mapq_scores, sramm_scores = line.strip().split('\t')
        stats[read_id] = (mapq_scores, sramm_scores)
        counter += 1
        if counter >= num_reads:
            mapq_hist, mapped_hist, uniq_rep_hist, unmapped_hist = batch_parse_stats(pool, chunk_size, stats, mapq_max,
                                                                                     num_score_bins, mapq_hist,
                                                                                     mapped_hist, uniq_rep_hist,
                                                                                     unmapped_hist)
            stats.clear()
            counter = 0
    mapq_hist, mapped_hist, uniq_rep_hist, unmapped_hist = batch_parse_stats(pool, chunk_size, stats, mapq_max,
                                                                             num_score_bins, mapq_hist,
                                                                             mapped_hist, uniq_rep_hist,
                                                                             unmapped_hist)
    pool.close()
    pool.join()
    return mapq_hist, mapped_hist, uniq_rep_hist, unmapped_hist


def batch_parse_stats(pool, chunk_size, stats, mapq_max, num_score_bins, mapq_hist, mapped_hist, uniq_rep_hist,
                      unmapped_hist):
    """
    MP wrapper
    Manual batching and manual histogram generation (saves memory)
    """
    mapq_batch = []
    mapped_batch = []
    uniq_rep_batch = []
    unmapped_batch = []
    for r in pool.imap_unordered(parse_stats, [(k, v) for k, v in stats.items()], chunk_size):
        k, mapq, sramm = r
        for val in mapq:
            mapq_batch.append(int(val))
        mi, ur, um = sramm
        mapped_batch.append(mi)
        uniq_rep_batch.append(ur)
        unmapped_batch.append(um)
    temp = np.histogram(mapq_batch, bins=np.linspace(0, mapq_max, mapq_max + 1))
    mapq_hist += temp[0]
    batches = [np.histogram(batch, bins=np.linspace(0, 1, num_score_bins + 1)) for batch in
               [mapped_batch, uniq_rep_batch, unmapped_batch]]
    mapped_hist += batches[0][0]
    uniq_rep_hist += batches[1][0]
    unmapped_hist += batches[2][0]
    return mapq_hist, mapped_hist, uniq_rep_hist, unmapped_hist


def parse_stats(tup):
    """
    Parses read stats using MP after reading
    :param tup: (key, (aln_scores, mapq_scores, sramm_scores)) tuple
    :return: read_id, mapq value(s), sramm metric scores
    """
    # data[k] = [(int(x[0]), int(x[1])) for x in
    #                  [t.replace('(', '').replace(')', '').split(', ') for t in val[0].split('),(')]]
    k, val = tup
    mapq = [int(v) for v in val[0].split(',')]
    sramm = [float(v) for v in val[1].split(',')]
    return k, mapq, sramm


def generate_graphs(graph_data, mapq_max, num_score_bins):
    """
    Wrapper for generating graphs and writing graph data file
    Persistent manual histogram version
    :param graph_data: manual histograms from stat file
    :param mapq_max: max value of mapq for graph purposes
    :param num_score_bins: number of bins for the 3 SRAMM metrics.
    :return:
    """
    mapq_hist, mapped_hist, uniq_rep_hist, unmapped_hist = graph_data
    ofh = open(output_prefix + "graphs_hist_data.txt", 'w')
    ofh.write("mapq\t" + ",".join([str(int(v)) for v in mapq_hist]) + "\n")
    ofh.write("bin edges\t" + ",".join([str(np.around(v, 4)) for v in np.linspace(0, 1, num_score_bins + 1)]) + "\n")
    ofh.write("mapped_identity\t" + ",".join([str(int(v)) for v in mapped_hist]) + "\n")
    ofh.write("unique_repeat\t" + ",".join([str(int(v)) for v in uniq_rep_hist]) + "\n")
    ofh.write("unmapped\t" + ",".join([str(int(v)) for v in unmapped_hist]))
    ofh.close()
    graph_mapq(mapq_hist, mapq_max)
    graph_score(mapped_hist, num_score_bins, 'Mapped Identity Score Distribution', 'Mapped Identity Score',
                'Counts (log10)', output_prefix + 'graphs_mapped_score.png')
    graph_score(uniq_rep_hist, num_score_bins, 'Unique Repeat Score Distribution', 'Unique-Repeat Score',
                'Counts (log10)', output_prefix + 'graphs_uniq_rep_score.png')
    graph_score(unmapped_hist, num_score_bins, 'Unmapped Score Distribution', 'Unmapped Score',
                'Counts (log10)', output_prefix + 'graphs_unmapped_score.png')
    graph_subplots(mapq_hist, mapq_max, mapped_hist, uniq_rep_hist, unmapped_hist, num_score_bins)
    return


def graph_subplots(mapq_hist, mapq_max, mapped_hist, uniq_rep_hist, unmapped_hist, num_score_bins):
    """

    """
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import numpy as np
    fig = plt.figure(figsize=(18, 18))
    # fig.subplots_adjust(top=0.3, bottom=0.1, left=0.1, right=0.3, wspace=0.25, hspace=0.25)
    # MAPQ
    fig = plt.subplot(2, 2, 1)
    title = 'MAPQ Distribution'
    xlabel = 'MAPQ'
    ylabel = 'Counts (log10)'
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.hist([i for i in range(mapq_max)], bins=np.linspace(0, mapq_max, mapq_max + 1), weights=mapq_hist, rwidth=0.8)
    plt.yscale('log', nonpositive='clip')
    # MI
    fig = plt.subplot(2, 2, 2)
    plt.title('Mapped Identity Score Distribution')
    plt.xticks(np.linspace(0, 1, 11))
    plt.xlabel('Mapped Identity Score')
    plt.ylabel('Counts (log10)')
    plt.hist([i / num_score_bins + 1 / (num_score_bins * 2) for i in range(num_score_bins)], bins=np.linspace(0, 1, num_score_bins + 1),
             weights=mapped_hist, rwidth=0.7)
    plt.yscale('log', nonpositive='clip')
    # UR
    fig = plt.subplot(2, 2, 3)
    plt.title('Unique Repeat Score Distribution')
    plt.xticks(np.linspace(0, 1, 11))
    plt.xlabel('Unique-Repeat Score')
    plt.ylabel('Counts (log10)')
    plt.hist([i / num_score_bins + 1 / (num_score_bins * 2) for i in range(num_score_bins)], bins=np.linspace(0, 1, num_score_bins + 1),
             weights=uniq_rep_hist, rwidth=0.7)
    plt.yscale('log', nonpositive='clip')
    # UM
    fig = plt.subplot(2, 2, 4)
    plt.title('Unmapped Score Distribution')
    plt.xticks(np.linspace(0, 1, 11))
    plt.xlabel('Unmapped Score')
    plt.ylabel('Counts (log10)')
    plt.hist([i / num_score_bins + 1 / (num_score_bins * 2) for i in range(num_score_bins)], bins=np.linspace(0, 1, num_score_bins + 1),
             weights=unmapped_hist, rwidth=0.7)
    plt.yscale('log', nonpositive='clip')
    plt.savefig(output_prefix + 'graphs_subplots.png', dpi=600)
    plt.close()
    return


def graph_mapq(scores, mapq_max):
    """
    Graph MAPQ
    Uses 'agg' for non GUI or X sessions ie can generate and save graph in clusters
    :param scores: manual histogram
    :param mapq_max: max value of mapq for graph purposes
    :return:
    """
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import numpy as np
    title = 'MAPQ Distribution'
    xlabel = 'MAPQ'
    ylabel = 'Counts (log10)'
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.hist([i for i in range(mapq_max)], bins=np.linspace(0, mapq_max, mapq_max + 1), weights=scores, rwidth=0.8)
    plt.yscale('log', nonpositive='clip')
    plt.savefig(output_prefix + 'graphs_mapq.png', dpi=600)
    plt.close()
    return


def graph_score(scores, num_bins, title, xlabel, ylabel, output_file):
    """
    Graphing template for SRAMM stats.
    Uses 'agg' for non GUI or X sessions ie can generate and save graph in clusters
    :param scores: manual histogram
    :param num_bins: number of bins for the manual histogram
    :param title: title
    :param xlabel: xlabel
    :param ylabel: ylabel
    :param output_file: saved graph file
    :return:
    """
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.title(title)
    plt.xticks(np.linspace(0, 1, 11))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.hist([i / num_bins + 1 / (num_bins * 2) for i in range(num_bins)], bins=np.linspace(0, 1, num_bins + 1),
             weights=scores, rwidth=0.7)
    plt.yscale('log', nonpositive='clip')
    plt.savefig(output_file, dpi=600)
    plt.close()
    return


def run_stats(input_file, output_file, end, num_reads, num_cpu, ftime):
    """
    Wrapper for running stats computation
    """
    print('-STATS: ', np.around(time.time() - ftime, 2))
    print('--Starting parser: ', np.around(time.time() - ftime, 2))
    ofh = open(output_file, 'w')
    extension = input_file.split(".")[-1].lower()
    if extension == "bam":
        parse_alignment_file(input_file, ofh, end, num_reads, num_cpu)
    elif extension == "vg":
        parse_vg_file(input_file, ofh, end, num_reads, num_cpu)
    else:
        raise Exception('alignment file must be .bam for BAM format or .vg for VG refpos-table format')
    ofh.close()
    print('--Finished calculating stats: ', np.around(time.time() - ftime, 2))
    return


def run_filter(input_file, output_file, num_reads, num_cpu, score_ranges, ftime):
    """
    Wrapper for running filtering
    """
    print('-FILTER: ', np.around(time.time() - ftime, 2))
    print('--Filtering stats: ', np.around(time.time() - ftime, 2))
    if output_file == input_file:
        output_file = output_file.split(".txt")[0] + "_filt.txt"
        print("--Input and output are the same, new output file: {}".format(output_prefix + output_file))
    ofh = open(output_prefix + output_file, 'w')
    filter_stats_wrapper(num_reads, num_cpu, score_ranges, input_file, ofh)
    ofh.close()
    print('--Finished filtering stats: ', np.around(time.time() - ftime, 2))
    return output_prefix + output_file


def run_graphs(input_file, num_reads, num_cpu, mapq_max, num_score_bins, ftime):
    """
    Wrapper for running graphs generation
    """
    print('-GRAPHS: ', np.around(time.time() - ftime, 2))
    print('--Reading stats file: ', np.around(time.time() - ftime, 2))
    # read input file which has to be SRAMM results files
    graph_data = read_stat_file(input_file, num_reads, num_cpu, mapq_max, num_score_bins)
    print('--Generating graphs: ', np.around(time.time() - ftime, 2))
    generate_graphs(graph_data, mapq_max, num_score_bins)
    print('--Finished generating graphs: ', np.around(time.time() - ftime, 2))
    return


def main(args):
    global output_prefix
    # Process input arguments
    process_type = args.process_type
    input_file = args.input_file
    output_file = args.output_file
    output_prefix = args.output_prefix
    end = args.end
    num_reads = int(args.num_reads)
    num_cpu = int(args.num_cpu)
    map_range = [float(v) for v in args.map_range.split(',')]
    uniq_rep_range = [float(v) for v in args.uniq_rep_range.split(',')]
    unmap_range = [float(v) for v in args.unmap_range.split(',')]
    num_alns_range = [int(float(v)) for v in args.num_alns_range.split(',')]
    score_ranges = (map_range, uniq_rep_range, unmap_range, num_alns_range)
    mapq_max = int(args.mapq_max)
    num_score_bins = int(args.score_num_bins)
    ftime = time.time()
    # Check process type and run pipelines
    if process_type not in ['stats', 'filter', 'graphs']:
        raise Exception('process_type must be stats, filter, graph')
    # Stat pipeline
    elif process_type == 'stats':
        run_stats(input_file, output_file, end, num_reads, num_cpu, ftime)
        if args.filter_flag == "True":
            output_file = run_filter(output_file, output_file, num_reads, num_cpu, score_ranges, ftime)
        if args.graphs_flag == "True":
            run_graphs(output_file, num_reads, num_cpu, mapq_max, num_score_bins, ftime)
    # Filter pipeline
    elif process_type == 'filter':
        output_file = run_filter(input_file, output_file, num_reads, num_cpu, score_ranges, ftime)
        if args.graphs_flag == "True":
            run_graphs(output_file, num_reads, num_cpu, mapq_max, num_score_bins, ftime)
    # Graph pipeline
    elif process_type == 'graphs':
        run_graphs(input_file, num_reads, num_cpu, mapq_max, num_score_bins, ftime)
    print('-Finished: ', np.around(time.time() - ftime, 2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='SRAMM.py',
                                     usage='%(prog)s [options] process_type input_file',
                                     description='Short Read Alignment Mapping Metrics',
                                     epilog='Visit https://github.com/achon/sramm for more information.')
    parser.add_argument("process_type", help="Process type: stats, filter, graphs")
    parser.add_argument("input_file", help="Input file, depends on process and input_type")
    parser.add_argument("-output_file", help="Output stat file, default=sramm_output.txt", default="sramm_output.txt")
    parser.add_argument("-output_prefix", help="Output prefix for all files, default=''", type=str, default='')
    parser.add_argument("-end", help="single or paired end, default=paired", default="paired")
    parser.add_argument("-num_reads", help="Number of reads worth of alignments to keep in memory, default=10000",
                        default=10000)
    parser.add_argument("-num_cpu", help="Number of cpu, default=2", default=2)
    parser.add_argument("-graphs_flag", help="Generate graphs flag, default=False", default="False")

    parser.add_argument("-filter_flag", help="Filter flag, default=False", default="False")
    parser.add_argument("-mapq_max", help="Max MAPQ value in dataset, default=60", default=60)
    parser.add_argument("-score_num_bins", help="Number of bins for SRAMM scores for graphing", default=40)
    parser.add_argument("-map_range", help="Range of the map score to filter (accept), default=0.0,1.0",
                        default="0.0,1.0")
    parser.add_argument("-uniq_rep_range", help="Range of the uniq_rep score to filter (accept), default=0.0,1.0",
                        default="0.0,1.0")
    parser.add_argument("-unmap_range", help="Range of the unmap score to filter (accept), default=0.0,1.0",
                        default="0.0,1.0")
    parser.add_argument("-num_alns_range", help="Range of the number of alignments to filter (accept), default=0.0,100.0",
                        default="0.0,100.0")
    args = parser.parse_args()
    print(args)
    main(args)

