#!/usr/bin/python3

import argparse
import os
from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from multiprocessing import cpu_count
from subprocess import call


def get_gene():
    """
    You can edit gene.list if you want more or less chloroplast coding
    gene or other fragment as long as it was described in genbank file.
    Also you can directly add mitochrondria gene name after the list. Ensure
    you use correct genbank file."""
    wanted_gene_list = list()
    with open('gene.list', 'r') as raw:
        for line in raw.readlines():
            if line.startswith('#'):
                continue
            line = line.strip()
            wanted_gene_list.append(line)

    fragment = list()
    genomes = SeqIO.parse(sys.argv[1], 'gb')
    for genome in genomes:
        for feature in genome.features:
            if feature.type != 'gene' or 'gene' not in feature.qualifiers:
                continue
            position = list()
            if feature.location_operator != 'join':
                position.append([
                    int(feature.location.start),
                    int(feature.location.end)
                ])
            else:
                for i in feature.sub_features:
                    position.append([
                        int(i.location.start),
                        int(i.location.end)
                    ])
            for n, frag in enumerate(position):
                name = str(feature.qualifiers['gene'][0]).replace(' ', '_')
                if name not in wanted_gene_list:
                    continue
                sequence = str(genome.seq[frag[0]:frag[1]])
                if n > 0:
                    name = '{0}-{1}'.format(name, n+1)
                fragment.append([name, sequence])
    return fragment


def generate_query(fragment):
    """
    Generate fragment.fasta to BLAST. You can delete it freely."""
    handle = open('fragment.fasta', 'w')
    for gene in fragment:
        handle.write('>{0}\n{1}\n'.format(gene[0], gene[1]))
    handle.close()
    return 'fragment.fasta'

    call('makeblastdb -in {0} -out {1} -dbtype nucl'.format(contig_file,
                                                            contig_file),
         shell=True)

def blast(query_file, db_file, result='out/BLASTResult.xml'):
    """Here we use "max_hsps" to restrict only the first hsp, use
    "max_target_seqs" to restrict only the first matched sequence.
    """
    blast_result_file = os.path.join(arg.output, 'BlastResult.xml')
    cmd = nb(num_threads=cpu_count(),
             query=query_file,
             db=db_file,
             task='blastn',
             evalue=arg.evalue,
             max_hsps=1,
             max_target_seqs=1,
             outfmt=5,
             out=result)
    stdout, stderr = cmd()
    return blast_result_file


def parse(blast_result_file):
    parse_result = list()
    blast_result = SearchIO.parse(blast_result_file, 'blast-xml')
    for record in blast_result:
        if len(record) == 0:
            continue
        for i in record:
            parse_result.append([i[0][0].hit, i[0][0].query.id])
    return parse_result


def output(parse_result, contig_file, mode):
    contigs = SeqIO.parse(contig_file, 'fasta')
    annotated_contig = contig_file.split(sep='.')[0]
    handle = open('out/{0}_filtered.fasta'.format(annotated_contig), 'w')
    parse_result_d = {i[0].id: [] for i in parse_result}
    for record in parse_result:
        parse_result_d[record[0].id].append([record[0].seq, record[1]])
    for contig in contigs:
        if contig.id not in parse_result_d:
            continue
        if mode == '1':
            gene = parse_result_d[contig.id]
            for match in gene:
                new_seq = SeqRecord(
                    id='{0}|{1}|{2}'.format(
                        sys.argv[2].replace('.fasta', ''), 
                        match[1], 
                        contig.id),
                    description='',
                    seq=match[0]
                )
                gene_file = 'out/{0}-{1}.fasta'.format(annotated_contig,
                                                       match[1])
                handle_gene = open(gene_file, 'a')
                SeqIO.write(new_seq, handle_gene, 'fasta')
        else:
            SeqIO.write(contig, handle, 'fasta')
    handle.close()


def filter(contig_file, minium_length):
    contig_raw = SeqIO.parse(contig_file, 'fasta')
    contig_long = list()
    for contig in contig_raw:
        if(len(contig.seq) < minium_length):
            pass
        else:
            contig_long.append(contig)
    contig_long_file = '{0}-long.fasta'.format(contig_file)
    SeqIO.write(contig_long, contig_long_file, 'fasta')
    return contig_long_file


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
    >>>python3 annotate_contig.py reference_file contig_file mode

    Mode:
        1. Query contig against coding genes, then every contig will be
        annotated by gene name. You will only get fragment of contigs which
        was recognized via BLAST.
        matched in BLAST.
        2. Query contig in a whole genome. It only judge if contig was
        similiar to genome of given genbank file. In this mode, you get full
        length of contig.
        3. Query contigs in one file  against BLAST database generated from
        given reference fasta file. The most similiar sequence in contig will
        be write into files named as
        {reference_sequence_id}_{contig_id}.fasta.

    All results was set in 'out/'. Also you can set it by "-o".  """
    print(main.__doc__)
    if not os.path.exists('out'):
        os.makedirs('out')
    arg = argparse.ArgumentParser()
    arg.add_argument('-r', '--reference', dest='ref_file',
                     help='reference sequences file (fasta format)')
    arg.add_argument('-q', '--query', dest='query_file',
                     help='query file (fasta format)')
    arg.add_argument('-e', dest='evalue', default=1e-5, type=float,
                     help='evalue for BLAST')
    arg.add_argument('-m', '--mode', dest='mode', default='3',
                     help='query mode, see help info of program')
    arg.add_argument('-min_len', dest='minium_length', type=int,
                     default=300, help='minium length of contig')
    arg.add_argument('-o', '--output', dest='out', default='out',
                     help='output path')
    arg = arg.parse_args()
    global arg
    if arg.mode not in ('1', '2', '3'):
        raise ValueError('Bad command!\n')
    contig_file = filter(sys.argv[2], minium_length)
    if arg.mode == '1':
        fragment = get_gene()
        query_file = generate_query(fragment)
        xml_file = blast(query_file, contig_file)
        parse_result = parse(xml_file)
        output(parse_result, contig_file, mode)
    else:
        query_file = sys.argv[1].replace('.gb', '.fasta')
        SeqIO.convert(sys.argv[1], 'gb', query_file, 'fasta')
        xml_file = blast(query_file, contig_file)
        parse_result = parse(xml_file)
        output(parse_result, contig_file, arg.mode)


if __name__ == '__main__':
    main()
