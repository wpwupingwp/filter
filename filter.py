#!/usr/bin/python3
import argparse
import os
import re
from multiprocessing import cpu_count
from subprocess import run
from timeit import default_timer as timer

from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb


def safe(old):
    return re.sub(r'\W', '_', old)


def get_gene(ref_file):
    """
    You can edit gene.list if you want more or less chloroplast coding
    gene or other fragment as long as it was described in genbank file.
    Also you can directly add mitochrondria gene name after the list.
    Ensure you use correct genbank file."""
    wanted_gene_list = list()
    with open(args.gene_list, 'r') as raw:
        for line in raw.readlines():
            # omit commit line
            if line.startswith('#'):
                continue
            line = line.strip()
            wanted_gene_list.append(line)

    fragment = list()
    genomes = SeqIO.parse(ref_file, 'gb')
    for genome in genomes:
        for feature in genome.features:
            if feature.type != 'gene' or 'gene' not in feature.qualifiers:
                continue
            position = list()
            if feature.location_operator != 'join':
                position.append([int(feature.location.start),
                                 int(feature.location.end)])
            else:
                for i in feature.sub_features:
                    position.append([int(i.location.start),
                                     int(i.location.end)])
            for n, frag in enumerate(position):
                name = str(feature.qualifiers['gene'][0]).replace(' ', '_')
                if name not in wanted_gene_list:
                    continue
                sequence = str(genome.seq[frag[0]:frag[1]])
                if n > 0:
                    name = '{0}-{1}'.format(name, n+1)
                fragment.append([name, sequence])

    gene_file = os.path.join(args.out, 'fragment.fasta')
    with open(gene_file, 'w') as output_file:
        for gene in fragment:
            output_file.write('>{0}\n{1}\n'.format(gene[0], gene[1]))
    return gene_file


def blast(ref_file, query_file):
    """
    max_target_seqs was reported to having bug."""
    # MAX_TARGET_SEQS = 2
    db_file = os.path.join(args.out, ref_file)
    # hide output
    with open(os.path.join(args.out, 'log.txt'), 'w') as log:
        run('makeblastdb -in {0} -out {1} -dbtype nucl'.format(
            ref_file, db_file), stdout=log, shell=True)
    result = os.path.join(args.out, 'BlastResult.xml')
    cmd = nb(num_threads=cpu_count(),
             query=query_file,
             db=db_file,
             task='blastn',
             evalue=args.evalue,
             max_hsps=1,
             # max_target_seqs=MAX_TARGET_SEQS,
             outfmt=5,
             out=result)
    stdout, stderr = cmd()
    return result


def parse(blast_result_file):
    blast_result = SearchIO.parse(blast_result_file, 'blast-xml')
    querys = set()
    for query in blast_result:
        if len(query) == 0:
            continue
        if len(query) == 1:
            yield query[0][0]
            continue
        # Blast Result
        #   |_query
        #       |_hit
        #           |_hsp
        hits = list(query)
        hits.sort(key=lambda x: x[0].bitscore_raw, reverse=True)
        best_hsp = hits[0][0]
        second = hits[1][0]
        if best_hsp.bitscore_raw == second.bitscore_raw:
            if best_hsp.query_id in querys:
                pass
            else:
                querys.add(best_hsp.query_id)
                print(best_hsp.query_id, best_hsp.hit_id,
                      best_hsp.bitscore_raw, second.hit_id,
                      second.bitscore_raw)
                if args.same_out:
                    yield second, True
                    yield best_hsp, True
        else:
            # is_same
            yield best_hsp, False


def output(blast_result_file):
    filtered = os.path.join(
        args.out, os.path.splitext(args.query_file)[0]+'-filtered.fasta')
    handle = open(filtered, 'w')
    # {query_id+description: [hit_id, 0]}
    query_hit = dict()
    query_hit['miss'] = ['NOT_FOUND', 0]
    query_hit['bad'] = ['LOW_SCORE', 0]
    query_hit['all'] = ['TOTAL', 0]
    query_hit['same'] = ['SAME_SCORE', 0]
    for record, is_same in parse(blast_result_file):
        query_hit['all'][1] += 1
        if is_same:
            query_hit['same'][1] += 1
            with open(os.path.join(
                    args.out, args.query_file+'.same_score.fasta'),
                      'a') as same_score:
                SeqIO.write(record.query, same_score, 'fasta')
                continue
        if record.bitscore_raw < args.score:
            query_hit['bad'][1] += 1
            with open(os.path.join(
                    args.out, args.query_file+'.low_score.fasta'),
                      'a') as low_score:
                SeqIO.write(record.query, low_score, 'fasta')
            continue
        if record.query.description == '':
            query_hit[record.query.id] = [record.hit.id, 0]
        else:
            query_hit[record.query.id+' '+record.query.description] = [
                record.hit.id, 0]
    if args.fragment_out is not True:
        for record in SeqIO.parse(args.query_file, 'fasta'):
            # filter sequence missed in BLAST
            description = record.description
            if description in query_hit:
                query_hit[description][1] += 1
            # BLAST will remove ";" at the end of sequence id
            elif description[:-1] in query_hit:
                description = description[:-1]
                query_hit[description][1] += 1
            else:
                query_hit['miss'][1] += 1
                with open(os.path.join(
                        args.out, args.query_file+'.not_found.fasta'),
                          'a') as not_found:
                    SeqIO.write(record, not_found, 'fasta')
                continue

            info = '-'.join([query_hit[description][0],
                             os.path.splitext(args.query_file)[0]])
            output = safe(info)+'.fasta'
            record.id = info+'-'+record.description
            record.description = ''
            SeqIO.write(record, handle, 'fasta')
            # append rather overwrite
            with open(os.path.join(args.out,
                                   output), 'a') as output_file:
                SeqIO.write(record, output_file, 'fasta')
    else:
        for record in parse(blast_result_file):
            if record.bitscore_raw < args.score:
                continue
            if record.query.description == '':
                info = record.query.id
            else:
                info = record.query.id+' '+record.query.description
            query_hit[info][1] += 1
            info = '-'.join([record.hit.id+record.hit.description,
                             os.path.splitext(args.query_file)[0]])
            output = safe(info)+'.fasta'
            # record.query.id = ''
            record.query.description = info+'-'+record.query.description
            record.query.id = ''.join([record.query.description,
                                      record.query.id])
            record.query.description = ''

            # output to one file
            SeqIO.write(record.query, handle, 'fasta')
            with open(os.path.join(args.out,
                                   output), 'a') as output_file:
                # output seperately
                SeqIO.write(record.query, output_file, 'fasta')
    handle.close()
    statistics = filtered.replace('-filtered.fasta', '-count.csv')
    count = dict()
    for i in query_hit.values():
        try:
            count[i[0]] += i[1]
        except KeyError:
            count[i[0]] = i[1]
    with open(statistics, 'w') as stat:
        for line in count.items():
            stat.write('{0},{1}\n'.format(*line))
    print('#'*80)
    print('Seq_id\tnumber')
    for i in count.items():
        print('{}\t{}'.format(*i))


def main():
    arg = argparse.ArgumentParser()
    arg.add_argument('-r', dest='ref_file',
                     help='reference sequences file (fasta format)')
    arg.add_argument('-q', dest='query_file',
                     help='query file (fasta format)')
    arg.add_argument('-l', dest='gene_list', help='list of gene you want')
    arg.add_argument('-e', dest='evalue', default=1e-5,
                     type=float, help='evalue for BLAST')
    arg.add_argument('-o', dest='out', default='out', help='output path')
    arg.add_argument('-f', dest='fragment_out', action='store_true',
                     help='only output matched part of'
                     'query sequence rather than whole sequence')
    arg.add_argument('-s', dest='score', type=float, default=60.0,
                     help='BLAST score cutoff')
    arg.add_argument('-same_out', action='store_true',
                     help='output querys that have same score in various '
                     'reference')
    global args
    args = arg.parse_args()

    start = timer()
    if not os.path.exists(args.out):
        os.makedirs(args.out)
    if args.ref_file.endswith('.gb'):
        if args.gene_list is not None:
            args.ref_file = get_gene(args.ref_file)
        else:
            ref_file = os.path.splitext(args.ref_file)[0] + '.fasta'
            SeqIO.convert(args.ref_file, 'gb', ref_file, 'fasta')
            args.ref_file = ref_file
    xml_file = blast(args.ref_file, args.query_file)
    output(xml_file)
    # blast-xml is too big
    os.remove(xml_file)
    end = timer()
    print('='*80)
    print('Cost {:.3f} seconds.'.format(end-start))


if __name__ == '__main__':
    main()
