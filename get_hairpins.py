import parsers as ps
import handlers as hn
import funk as fk

gtf_parser = ps.Gtf('knowns_toy.gtf')
dro_genomes = ['dro' + x for x in ['Bia2', 'Ele2', 'Ere2', 'Eug2', 'Moj3', 'Per1', 'Rho2',
                                   'Sec1', 'Sim1', 'Suz1', 'Tak2', 'Vir3', 'Yak3']]
parser = ps.Maf('/Users/martin/Documents/drosophila/toy.maf', main_genome='dm6', genomes=dro_genomes)
lfolder = hn.RnaFold().Lfold(['-T', '25', '--noLP'])

print('getting trans_exons')
trans_exons = gtf_parser.get_trans_exon(annot='.')

print('retrieving relative start stops')
trans_rel_stst = fk.get_start_stops(trans_exons, gtf_parser)

print('reading maf')
parser.get_alignments()

print('extracting sequences')
trans_seqs = fk.get_seqs(trans_exons, parser)
parser = None

print('folding')
with open('hairpins', 'w') as fout:
    for trans in trans_seqs.keys():
        dm6_seq_gapped = trans_seqs[trans]['dm6']
        dm6_seq = dm6_seq_gapped.replace('-', '')
        hairpins = lfolder.compute(dm6_seq)
        for hairpin in hairpins:
            hairp_data = fk.analyze_hairpin(hairpin, dm6_seq_gapped, trans_seqs[trans])
            loc = fk.determine_location(hairpin, trans_rel_stst[trans])

            # check location!!!

            if hairp_data:
                fout.write('\t'.join([str(x) for x in [trans, hairpin['fold'], loc, len(hairpin['fold']), hairpin['energy']] + hairp_data]))
                fout.write('\n')