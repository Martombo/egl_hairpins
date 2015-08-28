import functions as fn
import handlers as hn
gp = fn.GappedSeq()

def consv_score(seqA, seqB):
    assert len(seqA) == len(seqB)
    score, leng = 0, 0
    for k in range(len(seqA)):
        if seqA[k] != '-':
            if seqA[k] == seqB[k]:
                score += 1
            leng += 1
    return float(score) / leng

def determine_location(hairpin, rel_stst_trans):
    mid_pos = hairpin['pos'] + (len(hairpin['fold']) / 2)
    if mid_pos < rel_stst_trans[0]:
        return '5UTR'
    elif hairpin['pos'] < rel_stst_trans[1]:
        return 'CDS'
    return '3UTR'

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

def check_hairpin(hairpin_fold, seq):
    """here are the necessary condition for a hairpin:
    is_simple
    20 <= len <= 80
    get_Aprimes() > 1
    """
    hairp_fold = fn.Fold(hairpin_fold, seq)
    simple = hairp_fold.is_simple()
    leng = 19 < hairp_fold.len < 81
    aPrimes = hairp_fold.get_Aprimes()
    return simple, leng, aPrimes

def perc_matches(fold):
    match_n = fold.count('(') + fold.count(')')
    return float(match_n) / len(fold)

def get_ortho_data(dro_seq, hairp_fold, dm6_hairp_gap):
    con_fold = hn.RnaFold().CostraintFold()
    seq, constr = gp.get_ortho_constr(dm6_hairp_gap, hairp_fold, dro_seq)
    if not seq:
        return '\t'.join([str(x) for x in [0, 0, 0, 0, False]])
    (fold, fold_score) = con_fold.compute(seq + '\n' + constr)
    check = check_hairpin(fold, seq)
    return '\t'.join([str(x) for x in [len(fold), fold_score, perc_matches(fold), len(check[2]), check[0]]])

def analyze_hairpin(hairpin, dm6_seq_gapped, trans_seqs_trans):
    (start, stop) = gp.get_gapped_i(dm6_seq_gapped, hairpin['pos']-1, len(hairpin['fold']))
    dm6_hairp_gap = dm6_seq_gapped[start:stop]
    dm6_hairp = dm6_hairp_gap.replace('-','')
    check = check_hairpin(hairpin['fold'], dm6_hairp)
    if check[0] and check[1] and len(check[2]) > 0:
        consv, ortho_data = [], ''
        for dro in sorted(trans_seqs_trans.keys()):
            if dro == 'dm6':
                continue
            dro_seq = trans_seqs_trans[dro]
            ortho_data += get_ortho_data(dro_seq[start:stop], hairpin['fold'], dm6_hairp_gap) + '\t'
            consv.append(consv_score(dm6_hairp_gap, dro_seq[start:stop]))
        cScore = sum(consv) / len(consv)
        return '\t'.join(str (x) for x in [dm6_hairp, len(check[2]), cScore, ortho_data])
