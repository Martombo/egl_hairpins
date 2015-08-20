import parsers as ps
import functions as fn

gtf_parser = ps.Gtf('/Users/martin/Dropbox/utils/drosophila/Drosophila_melanogaster.BDGP6.81toy.gtf')

trans_exons = gtf_parser.get_trans_exon(annot='.')
trans_exons_manager = fn.TransExons(trans_exons)
trans_start_stops = gtf_parser.get_start_stop_codons(annot='.')

trans_rel_start_stops = {}
for trans in trans_exons.keys():
    if trans in trans_start_stops:
        (start, stop) = trans_start_stops[trans]
        trans_rel_start_stops[trans] = trans_exons_manager.rel_pos_trans(trans, [start[0], stop[0]])

