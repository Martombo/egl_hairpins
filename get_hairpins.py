import parsers as ps
import functions as fn
import handlers as hn

gtf_parser = ps.Gtf('knowns_toy.gtf')
dro_genomes = ['dro' + x for x in ['Bia2', 'Ele2', 'Ere2', 'Eug2', 'Moj3', 'Per1', 'Rho2',
                                   'Sec1', 'Sim1', 'Suz1', 'Tak2', 'Vir3', 'Yak3']]
parser = ps.Maf('/Users/martin/Documents/drosophila/toy.maf', main_genome='dm6', genomes=dro_genomes)


trans_exons = gtf_parser.get_trans_exon(annot='.')
trans_exons_manager = fn.TransExons(trans_exons)
trans_start_stops = gtf_parser.get_start_stop_codons(annot='.')
lfolder = hn.RnaFold().Lfold()

trans_rel_stst = {}
for trans in trans_exons.keys():
    if trans in trans_start_stops:
        start_stop = trans_start_stops[trans]
        if len(start_stop) != 2:
            continue
        trans_rel_stst[trans] = trans_exons_manager.rel_pos_trans(trans, [start_stop[0][0], start_stop[1][0]])

parser.get_alignments()

trans_seqs = {}

for trans in trans_exons.keys():
    seqs = {}
    for exon in trans_exons[trans]:
        region = parser.get_region(exon[1], exon[2], exon[3])
        for key, obj in region.items():
            seqs = fn.add2dict(seqs, key, obj)
    for key, obj in seqs.items():
        seqs[key] = ''.join(obj)
    trans_seqs[trans] = seqs

parser = None

for trans in trans_seqs.keys():
    dm6_seq = trans_seqs[trans]['dm6']
    hairpins = lfolder.compute(dm6_seq.replace('-', ''), temp=25, noLP=True)
    for hairpin in hairpins:
        pos = hairpin['pos'] - 1
        seq = dm6_seq.replace('-','')[pos : pos + len(hairpin['fold'])]
        folder = fn.Fold(hairpin['fold'], seq)
        print(hairpin['fold'] + '\n' + seq)