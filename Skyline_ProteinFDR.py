from itertools import islice
import cmath
import math
from statsmodels.stats.multitest import multipletests

def extract_feature_info(infile):
    target_pros = {}
    decoy_peps = {}
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split(',')
            protein = split_i[8]
            pep = split_i[6]
            mod_pep = split_i[7]
            mprophet_score = split_i[10]
            qvalue = split_i[12]
            pvalue = split_i[11]
            rt = split_i[3]
            if protein != 'Decoys':
                if float(qvalue) < 0.01:
                    if protein not in target_pros:
                        target_pros[protein] = [pep + '@' + mod_pep + '@' + mprophet_score + '@' + pvalue + '@' + qvalue + '@' + rt]
                    else:
                        target_pros[protein].append(pep + '@' + mod_pep + '@' + mprophet_score + '@' + pvalue + '@' + qvalue + '@' + rt)

            elif protein == 'Decoys':
                if pep not in decoy_peps:
                    decoy_peps[pep] = [mprophet_score + '@' + pvalue + '@' + qvalue]
                else:
                    decoy_peps[pep].append(mprophet_score + '@' + pvalue + '@' + qvalue)

    return target_pros, decoy_peps

def get_peptide_vals(peps, pep_type):
    dicts = {}
    for peptide, vals in peps.items():
        for idx, scores in enumerate(vals):
            mproph_score = scores.split('@')[0]
            pval = scores.split('@')[1]
            qval = scores.split('@')[2]
            if pep_type == 'decoy':
                if 'decoy' in peptide:
                    if mproph_score not in dicts:
                        dicts[mproph_score] = [peptide]
                    else:
                        dicts[mproph_score].append(peptide)
                        
            elif pep_type == 'target':
                if mproph_score not in dicts:
                    dicts[mproph_score] = [peptide]
                else:
                    dicts[mproph_score].append(peptide)

    return dicts

def correct_pvals(pvals):

    corrected_vals = multipletests(pvals, alpha=0.01, method='fdr_tsbh', is_sorted=False, returnsorted=False)

    return corrected_vals

def calculate_protein_qval(target_pros, decoy_peps):            
    c = 0
    protein_qvals = {}
    for k, v in target_pros.items():
        peps = {}
        for p in v:
            pep = p.split('@')[0]
            decoy_pep = pep[0:len(pep)-1][::-1] + pep[-1]
            if decoy_pep in decoy_peps:
                for decoy_vals in decoy_peps[decoy_pep]:
                    if decoy_pep + '_decoy' not in peps:
                        peps[decoy_pep + '_decoy'] = [decoy_vals]
                    else:
                        peps[decoy_pep + '_decoy'].append(decoy_vals)
            
            if pep not in peps:
                peps[pep] = [p.split('@')[2] + '@' + p.split('@')[3] + '@' + p.split('@')[4]]
            else:
                peps[pep].append(p.split('@')[2] + '@' + p.split('@')[3] + '@' + p.split('@')[4])


        decoy_scores = get_peptide_vals(peps, 'decoy')
        target_scores = get_peptide_vals(peps, 'target')

        c += 1
        max_target_score = float(max(list(target_scores)))
        try:
            max_decoy_score = float(max(list(decoy_scores)))
        except:
            print (decoy_scores)
            
        competition_ratio = abs(max_target_score/max_decoy_score)
        adjusted_qvals = []
        for peptide, vals in peps.items():
            if 'decoy' not in peptide:
                for idx, scores in enumerate(vals):
                    mproph_score = scores.split('@')[0]
                    pval = float(scores.split('@')[1])
                    qval = float(scores.split('@')[2])
                    
                    #print (k, peptide, mproph_score, pval, qval, competition_ratio, competition_ratio*qval)

                    adjusted_qvals.append(competition_ratio*qval)

        #print (k, adjusted_qvals)
        if min(adjusted_qvals) not in protein_qvals:
            protein_qvals[min(adjusted_qvals)] = [k]
        else:
            protein_qvals[min(adjusted_qvals)].append(k)

    return protein_qvals, c


def protein_fdr(Skyline_mProphet_features, qval_correction):
    '''Extract the target proteins and their corresponding peptides with the mProphet Feature score, precursor q-value and p-value.'''
    target_proteins, decoy_peptides = extract_feature_info(Skyline_mProphet_features)

    '''The protein q-value will be calculated by calculating the competetion ration between highly scored target and decoy peptides of each protein. This competetion score will be multiplied with the peptide q-value to get the Protein q-value.'''
    protein_qvals, target_pro_count = calculate_protein_qval(target_proteins, decoy_peptides)

    pg = 0
    PGs = {}
    if qval_correction == True:
        ''' The protein q-value will be further subjected to FDR correction based on "fdr_tsbh" (two stage fdr correction (non-negative)) method.'''
        ''' This is under trial.'''
        b, qvals_corrected, alphacSidak, alphacBonf = correct_pvals(list(protein_qvals))
        print (f'INFO: Total number of q-value corrected proteins are {len(qvals_corrected)} as same as that of total number of proteins ({len(protein_qvals)}) with q-value.')
        
        for idx, qvals in enumerate(list(protein_qvals)):
            #print (idx, qvals, b[idx], qvals_corrected[idx], protein_qvals[qvals])
            if qvals_corrected[idx] < 0.01:
                for p in protein_qvals[qvals]:
                    pg += 1
                    pro_pep = {}
                    for t in target_proteins[p]:
                        if float(t.split('@')[4]) < 0.01:
                            if t.split('@')[1] not in pro_pep:
                                pro_pep[t.split('@')[1]] = [t]
                            else:
                                pro_pep[t.split('@')[1]].append(t)
                        
                    #print (p, max(scores), score[max(scores)])
                    for pep, _pep in pro_pep.items():
                        scores = {}
                        for _p in _pep:
                            if _p.split('@')[4] not in scores:
                                scores[_p.split('@')[4]] = [_p]
                            else:
                                scores[_p.split('@')[4]].append(_p)

                        #print (pep, max(scores), scores[max(scores)])
                        PGs[p + '@' + str(round(qvals_corrected[idx],4)) + '@' + scores[max(scores)][0].split('@')[1] + '@' + scores[max(scores)][0].split('@')[-1]] = 1
                        
    elif qval_correction == False:
        print (f'INFO: The q-value correction was not performed.')
        for idx, qvals in enumerate(list(protein_qvals)):
            if qvals < 0.01:
                for p in protein_qvals[qvals]:
                    pg += 1
                    pro_pep = {}
                    for t in target_proteins[p]:
                        if float(t.split('@')[4]) < 0.01:
                            if t.split('@')[1] not in pro_pep:
                                pro_pep[t.split('@')[1]] = [t]
                            else:
                                pro_pep[t.split('@')[1]].append(t)
                            
                    for pep, _pep in pro_pep.items():
                        scores = {}
                        for _p in _pep:
                            if _p.split('@')[4] not in scores:
                                scores[_p.split('@')[4]] = [_p]
                            else:
                                scores[_p.split('@')[4]].append(_p)

                        #print (pep, max(scores), scores[max(scores)])
                        PGs[p + '@' + str(round(qvals,4)) + '@' + scores[max(scores)][0].split('@')[1] + '@' + scores[max(scores)][0].split('@')[-1]] = 1
                
    pep_count = len(PGs)

    print (f"INFO: There are {pg} proteins with Protein q-value < 0.01 out of {target_pro_count} target proteins identified")
    print (f"INFO: There are {pep_count} peptides with FDR 1% cutoff supporting in the detection of {pg} high confidence proteins")

    output = [k.split('@') for k, v in PGs.items()]
    outfile = "{0}_ProteinGroups.txt".format(infile.rstrip('_mProphet_Features.csv'))
    with open(outfile, 'w') as outf:
        outf.write('Protein Accession\tProtein Q-value\tPeptide Sequence\tRT(min)\n')
        outf.writelines('\t'.join(i) + '\n' for i in output)

infile = "../M_tuberculosis_H37Rv/Mtb_DIA_Proteomics_SpecLib_Search_062023_mProphet_features.csv"

if __name__ == '__main__':
    protein_fdr(infile, False)

