#!/usr/bin/python3

import argparse
import os
from Bio import SearchIO, SeqIO, SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from multiprocessing import cpu_count
from subprocess import call
from tempfile import mkdtemp


def filter_length():
    long_contig = os.path.join(args.tmp, args.contig_file)
    with open(long_contig, 'w') as output_file:
        for contig in SeqIO.parse(args.contig_file, 'fasta'):
            if len(contig.seq) < args.min_length:
                pass
            else:
                SeqIO.write(contig, output_file, 'fasta')
    return long_contig


def get_gene(ref_file):
    """
    You can edit gene.list if you want more or less chloroplast coding
    gene or other fragment as long as it was described in genbank file.
    Also you can directly add mitochrondria gene name after the list. Ensure
    you use correct genbank file."""
    wanted_gene_list = list()
    with open(args.gene_list, 'r') as raw:
        for line in raw.readlines():
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
    Here it uses "max_hsps" to restrict only the first hsp,
    uses "max_target_seqs" to restrict only the first matched sequence."""
    db_file = os.path.join(args.tmp, ref_file)
    call('makeblastdb -in {0} -out {1} -dbtype nucl'.format(
        ref_file, db_file), shell=True)
    result = os.path.join(args.tmp, 'BlastResult.xml')
    cmd = nb(num_threads=cpu_count(),
             query=query_file,
             db=db_file,
             task='blastn',
             evalue=args.evalue,
             # to be continue
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
            # to be continued
            parse_result.append([i[0][0].hit, i[0][0].query])
    return parse_result


def output(parse_result):
    # to be continued
    annotated_contig = os.path.join(
        args.out, args.contig_file.replace('.fasta', '')+'_filter.fasta')
    handle = open(annotated_contig, 'w')
    if args.gene_list is not None:
        parse_result_d = {i[0].id: [] for i in parse_result}
        for record in parse_result:
            parse_result_d[record[0].id].append([record[0].seq, record[1]])

    contigs = SeqIO.parse(args.query_file, 'fasta')
    for contig in contigs:
        if contig.id not in parse_result_d:
            continue
        if args.fragment_out is not True:
            SeqIO.write(contig, handle, 'fasta')
        else:
            gene = parse_result_d[contig.id]
            for match in gene:
                new_seq = SeqRecord(
                    id='{0}|{1}|{2}'.format(
                        # to be continued
                        args.ref_file.split(sep='.')[0],
                        match[1],
                        contig.id),
                    description='',
                    seq=match[0]
                )
                gene_file = os.path.join(args.out, match[1]+'.fasta')
                # gene_file = 'out/{0}-{1}.fasta'.format(
                # annotated_contig, match[1])
                with open(gene_file, 'a') as gene_out:
                    SeqIO.write(new_seq, gene_out, 'fasta')
    handle.close()


def output2(parse_result):
    # to be continued
    filtered = os.path.join(args.out, args.query_file.replace(
        '.fasta', '')+'_filter.fasta')
    handle = open(filtered, 'w')
    if args.fragment_out is not True:
        # {query_id:hit_id}
        query_hit = {i[1].id: i[0].id for i in parse_result}
        for record in SeqIO.parse(args.query_file, 'fasta'):
            # filter sequence missed in BLAST
            if record.id not in query_hit:
                continue
            info = '-'.join([args.query_file.replace('.fasta', ''),
                             query_hit[record.id]])
            output = info+'.fasta'
            record.description = record.description+'-'+info
            SeqIO.write(record, handle, 'fasta')
            # append rather overwrite
            with open(os.path.join(args.out, output), 'a') as output_file:
                SeqIO.write(record, output_file, 'fasta')
    else:
        for record in parse_result:
            info = args.query_file.replace('.fasta', '')+'-'+record[0].id
            output = info+'.fasta'
            record[1].description = record[1].description+'-'+info
            SeqIO.write(record[1], handle, 'fasta')
            with open(os.path.join(args.out, output), 'w') as output_file:
                SeqIO.write(record[1], output_file, 'fasta')
    handle.close()


def main():
    """
    This program will annotate contigs from assembled sequences  according to
    given genbank file which describes a complete chloroplast genome or fasta
    format as reference sequence. The genbank file may contains one or more
    genomes.

    Edit gene.list if you want to annotate mitochrondria contigs.
    Notice that contig shorter than 300 bp will be ignored. You can change the
    minium length as you wish.

    Usage:
    >>>python3 annotate_contig.py -r reference_file -q contig_file

    Mode:
        1. Query contig against coding genes extract from given genbank file,
        then every contig will be annotated by gene name. You will only get
        fragment of contigs which was matched against reference sequences by
        BLAST. For instance :
        >>> python3 annotated_contig.py -r ref.gb -q contigs.fasta -l gene.list

        2. Query contig in  whole sequence. It only judges if contig was
        similiar to given reference sequence. In this mode, you get full
        length of contig.
        3. Query contigs in one file  against BLAST database generated from
        given reference fasta file. The most similiar sequence in contig will
        be write into files named as
        {contig_filename}_{reference_sequence_id}.fasta.
        """
    print(main.__doc__)
    arg = argparse.ArgumentParser()
    arg.add_argument('-r', dest='ref_file',
                     help='reference sequences file (fasta format)')
    arg.add_argument('-q', dest='query_file',
                     help='query file (fasta format)')
    arg.add_argument('-l', dest='gene_list', help='list of gene you want')
    arg.add_argument('-e', dest='evalue', default=1e-5,
                     type=float, help='evalue for BLAST')
    arg.add_argument('mode', type=int, choices=(1, 2, 3),
                     help='query mode, see help info of program')
    arg.add_argument('-min_len', dest='min_length', type=int,
                     default=10, help='minium length of contig')
    arg.add_argument('-o', dest='out', default='out',
                     help='output path')
    arg.add_argument('-f', dest='fragment_out', action='store_true',
                     help='only output matched part of'
                     'query sequence rather than whole sequence')
    arg.add_argument('-tmpdir', dest='tmp', default=mkdtemp(),
                     help='temporary directory')
    global args
    args = arg.parse_args()

    if not os.path.exists(args.out):
        os.makedirs(args.out)
    # filter length
    try:
        args.query_file = filter_length()
    except:
        arg.print_help()
    if args.ref_file.endswith('.gb'):
        if args.gene_list is not None:
            args.ref_file = get_gene(args.query_file)
        else:
            query_file = args.ref_file.replace('.gb', '.fasta')
            SeqIO.convert(args.ref_file, 'gb', query_file, 'fasta')
            args.query_file = query_file
    xml_file = blast(args.ref_file, args.query_file)
    parse_result = parse(xml_file)
    output2(parse_result)


if __name__ == '__main__':
    main()
