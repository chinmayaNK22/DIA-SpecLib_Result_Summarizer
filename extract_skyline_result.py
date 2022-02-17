from itertools import islice
import argparse

#infile = "NMO_DIA_25Da_PI_6_3_0.02Da_20ppm_spec_lib_search_11042021_Result.tsv"

parser = argparse.ArgumentParser(description='''Extract and summarize Skyline results from DIA data analysis against spectral library''')

parser.add_argument('infile', metavar='-ip', type=str, nargs='+', help='Skyline output from DIA-Spectral Library search')

args = parser.parse_args()

def get_header_idx(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            split_i = i.rstrip().split('\t')
            pep_idx = split_i.index("Peptide")
            pro_idx = split_i.index("Protein Name")
            qvalue_idx = split_i.index("Detection Q Value")
            z = split_i.index("Precursor Charge")
            raw_file = split_i.index("Replicate")
            miss_cleave_idx = split_i.index("Missed Cleavages")
            mz = split_i.index("Precursor Mz")
            rt = split_i.index("Peptide Retention Time")
            return pep_idx, pro_idx, qvalue_idx, z, miss_cleave_idx, mz, rt, raw_file

def ext_skyline_results(infile):            
    dicts = {}
    dia_file = {}
    a = get_header_idx(infile)
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            if split_i[a[0]] + '@' + split_i[a[1]] + '@' + split_i[a[3]] not in dicts:
                dicts[split_i[a[0]] + '@' + split_i[a[1]] + '@' + split_i[a[3]]] = [split_i]
            else:
                dicts[split_i[a[0]] + '@' + split_i[a[1]] + '@' + split_i[a[3]]].append(split_i)

            if split_i[a[-1]] not in dia_file:
                dia_file[split_i[a[-1]]] = [split_i]
            else:
                dia_file[split_i[a[-1]]].append(split_i)

    pep = {k.split('@')[0]:k.split('@')[0] + '@' + k.split('@')[2] for k, v in dicts.items() for j in v if j[3] != "#N/A" if float(j[a[2]]) < 0.01}
    pro = {k.split('@')[1]:k.split('@')[0] + '@' + k.split('@')[2] for k, v in dicts.items() for j in v if j[3] != "#N/A" if float(j[a[2]]) < 0.01}
    summary = [str(len(dicts)), str(len(pep)), str(len(pro))]
    summary_head = ['No. of Precursors','No. of Peptides','No. of Proteins']
    output = []
    z = {}
    mmc = {}
    output1 = []
    for k, v in dicts.items():
        charge = ''.join(m for m, n in {j[a[3]]:j for j in v}.items())
        miss_cleave = ''.join(m for m, n in {j[a[4]]:j for j in v}.items())
        fdr = min(j[a[2]] for j in v)
        if charge not in z:
            z[charge] = [k]
        else:
            z[charge].append(k)

        if miss_cleave not in mmc:
            mmc[miss_cleave] = [k]
        else:
            mmc[miss_cleave].append(k)
        output.append(k.split('@') + [charge] + [miss_cleave] + [fdr])

    raw_file = []
    for k, v in dia_file.items():
        raw_file.append(k)
        precursors = {}
        for j in v:
            if j[a[0]] + '@' + j[a[1]] + '@' + j[a[5]] not in precursors:
                precursors[j[a[0]] + '@' + j[a[1]] + '@' + j[a[5]] + '@' + j[a[6]]] = [j]
            else:
                precursors[j[a[0]] + '@' + j[a[1]] + '@' + j[a[5]] + '@' + j[a[6]]].append(j)
                
        for m, n in precursors.items():
            fdr = ';'.join(j[a[2]] for j in n)
            output1.append([k] + m.split('@') + [fdr])

    write1 = open("{0}_Charge_state_summary.txt".format(infile.rstrip('.tsv')), 'w')
    for k, v in z.items():
        write1.write(k + '\t' + str(len(v)) + '\n')
    write1.close()

    write2 = open("{0}_missed_clavage_summary.txt".format(infile.rstrip('.tsv')), 'w')
    for k, v in mmc.items():
        write2.write(k + '\t' + str(len(v)) + '\n')
    write2.close()
                  
    outfile = "{0}_peptides_and_proteins.txt".format(infile.rstrip('.tsv'))
    with open(outfile, 'w') as outf:
        outf.write('Protein\tPeptide\tPrecursor\tCharge\tMissed cleave\tQ-value\n')
        outf.writelines('\t'.join(i) + '\n' for i in output)

    outfile1 = "{0}_raw_file_specific_peptides_and_proteins.txt".format(infile.rstrip('.tsv'))
    with open(outfile1, 'w') as outf1:
        outf1.write('Raw File\tProtein\tPeptide\tPrecursor\tRT\tQ-values\n')
        outf1.writelines('\t'.join(i) + '\n' for i in output1)

    summrayfile = "{0}_summary.txt".format(infile.rstrip('.tsv'))
    with open(summrayfile, 'w') as sumf:
        sumf.write('\t'.join(summary_head) + '\n')
        sumf.write('\t'.join(summary) + '\n')

if __name__== "__main__":
    ext_skyline_results(args.infile[0])
