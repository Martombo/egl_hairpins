import functions as fn
import handlers as hn
gp = fn.GappedSeq()

def determine_location(hairpin, rel_stst_trans):
    mid_pos = hairpin['pos'] + (len(hairpin['fold']) / 2)
    if mid_pos < rel_stst_trans[0]:
        return '5\'UTR'
    elif hairpin['pos'] < rel_stst_trans[1]:
        return 'CDS'
    return '3\'UTR'

def get_start_stops(trans_exons, gtf_parser):
    trans_exons_manager = fn.TransExons(trans_exons)
    trans_start_stops = gtf_parser.get_start_stop_codons(annot='.')
    trans_rel_stst = {}
    for trans in trans_exons.keys():
        if trans in trans_start_stops:
            start_stop = trans_start_stops[trans]
            if len(start_stop) != 2:
                continue
            trans_rel_stst[trans] = trans_exons_manager.rel_pos_trans(trans, [start_stop[0][0], start_stop[1][0]])
    return trans_rel_stst

def get_seqs(trans_exons, parser):
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
    return trans_seqs

def is_hairpin_valid(hairpin_fold, seq):
    """here are the necessary condition for a hairpin:
    20 <= len <= 80
    is_simple
    get_Aprimes() > 1
    """
    hairp_fold = fn.Fold(hairpin_fold, seq)
    if not hairp_fold.is_simple():
        return False
    if hairp_fold.len < 20 or hairp_fold.len > 80:
        return False
    aPrimes = hairp_fold.get_Aprimes()
    if len(aPrimes) > 1:
        return aPrimes

def perc_matches(fold):
    match_n = fold.count('(') + fold.count(')')
    return float(match_n) / len(fold)

def get_ortho_data(dro_seq, hairp_fold, dm6_hairp_gap):
    con_fold = hn.RnaFold().CostraintFold()
    seq, constr = gp.get_ortho_constr(dm6_hairp_gap, hairp_fold, dro_seq)
    if not seq:
        return
    (fold, fold_score) = con_fold.compute(seq + '\n' + constr)
    aPrimes = is_hairpin_valid(fold, seq)

def analyze_hairpin(hairpin, dm6_seq_gapped, trans_seqs_trans):
    (start, stop) = gp.get_gapped_i(dm6_seq_gapped, hairpin['pos']-1, len(hairpin['fold']))
    dm6_hairp_gap = dm6_seq_gapped[start:stop]
    dm6_hairp = dm6_hairp_gap.replace('-','')
    aPrimes = is_hairpin_valid(hairpin['fold'], dm6_hairp)
    if aPrimes:
        ortho_data = {}
        for dro in trans_seqs_trans:
            dro_seq = trans_seqs_trans[dro]
            ortho_data[dro] = get_ortho_data(dro_seq[start:stop], hairpin['fold'], dm6_hairp_gap)
        return [dm6_hairp, perc_matches(hairpin['fold']), ortho_data]
    return
