#!/usr/bin/python3
import argparse
import os
from Bio import SearchIO, SeqIO, SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from multiprocessing import cpu_count
from subprocess import call
from tempfile import mkdtemp


def filter_length():
    long_contig = os.path.join(args.tmp, args.query_file)
    with open(long_contig, 'w') as output_file:
        for contig in SeqIO.parse(args.query_file, 'fasta'):
            if len(contig.seq) < args.min_length:
                pass
            else:
                SeqIO.write(contig, output_file, 'fasta')
    return long_contig


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

    gene_file = os.path.join(args.tmp, 'fragment.fasta')
    with open(gene_file, 'w') as output_file:
        for gene in fragment:
            output_file.write('>{0}\n{1}\n'.format(gene[0], gene[1]))
    return gene_file


def blast(ref_file, query_file):
    """
    Here it uses "max_hsps" to restrict only the first hsp, uses
    "max_target_seqs" to restrict only the first matched sequence."""
    db_file = os.path.join(args.tmp, ref_file)
    call('makeblastdb -in {0} -out {1} -dbtype nucl'.format(
        ref_file, db_file), shell=True)
    result = os.path.join(args.tmp, 'BlastResult.xml')
    cmd = nb(num_threads=cpu_count(),
             query=query_file,
             db=db_file,
             task='blastn',
             evalue=args.evalue,
             max_hsps=1,
             max_target_seqs=1,
             outfmt=5,
             out=result)
    stdout, stderr = cmd()
    return result


def parse(blast_result_file):
    parse_result = list()
    blast_result = SearchIO.parse(blast_result_file, 'blast-xml')
    for record in blast_result:
        if len(record) == 0:
            continue
        for i in record:
            parse_result.append([i[0][0].hit, i[0][0].query])
    return parse_result


def output(parse_result):
    filtered = os.path.join(args.out, os.path.basename(
        #  /tmp/tmpabcdef/out.fasta -> out.fasta to avoid wrong output path
        args.query_file.replace('.fasta', '')+'_filter.fasta'))
    handle = open(filtered, 'w')
    if args.fragment_out is not True:
        # {query_id:hit_id}
        query_hit = {i[1].id: i[0].id for i in parse_result}
        for record in SeqIO.parse(args.query_file, 'fasta'):
            # filter sequence missed in BLAST
            if record.id not in query_hit:
                continue
            info = '-'.join([os.path.basename(
                args.query_file.replace('.fasta', '')), query_hit[record.id]])
            output = info+'.fasta'
            record.description = record.description+'-'+info
            SeqIO.write(record, handle, 'fasta')
            # append rather overwrite
            with open(os.path.join(args.out,
                                   output), 'a') as output_file:
                SeqIO.write(record, output_file, 'fasta')
    else:
        for record in parse_result:
            info = os.path.basename(
                args.query_file.replace('.fasta', '')+'-'+record[0].id)
            output = info+'.fasta'
            record[1].description = record[1].description+'-'+info
            # output to one file
            SeqIO.write(record[1], handle, 'fasta')
            with open(os.path.join(args.out,
                                   output), 'w') as output_file:
                # output seperately
                SeqIO.write(record[1], output_file, 'fasta')
    handle.close()


def main():
    arg = argparse.ArgumentParser()
    arg.add_argument('-r', dest='ref_file',
                     help='reference sequences file (fasta format)')
    arg.add_argument('-q', dest='query_file',
                     help='query file (fasta format)')
    arg.add_argument('-l', dest='gene_list', help='list of gene you want')
    arg.add_argument('-e', dest='evalue', default=1e-5,
                     type=float, help='evalue for BLAST')
    arg.add_argument('-min_len', dest='min_length', type=int,
                     default=10, help='minium length of contig')
    arg.add_argument('-o', dest='out', default='out', help='output path')
    arg.add_argument('-f', dest='fragment_out', action='store_true',
                     help='only output matched part of'
                     'query sequence rather than whole sequence')
    arg.add_argument('-tmpdir', dest='tmp', default=mkdtemp(),
                     help='temporary directory')
    global args
    args = arg.parse_args()

    if not os.path.exists(args.out):
        os.makedirs(args.out)
    if args.ref_file.endswith('.gb'):
        if args.gene_list is not None:
            args.ref_file = get_gene(args.ref_file)
        else:
            ref_file = args.ref_file.replace('.gb', '.fasta')
            SeqIO.convert(args.ref_file, 'gb', ref_file, 'fasta')
            args.ref_file = ref_file
    # filter length
    args.query_file = filter_length()
    xml_file = blast(args.ref_file, args.query_file)
    parse_result = parse(xml_file)
    output(parse_result)
    print('Done.')


if __name__ == '__main__':
    main()
