from itertools import islice
import argparse

parser = argparse.ArgumentParser(description='''Extract and summarize Skyline results from DIA data analysis against spectral library''')

parser.add_argument('infile', metavar='-ip', type=str, nargs='+', help='Skyline output from DIA-Spectral Library search in tab delimitted format (txt/tsv)')

args = parser.parse_args()

def get_header_idx(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            split_i = i.rstrip().split('\t')
            pep_idx = split_i.index("Peptide")
            try:
                mod_pep_idx = split_i.index("Modified Sequence")
            except:
                mod_pep_idx = split_i.index("Peptide")
            pro_idx = split_i.index("Protein Name")
            qvalue_idx = split_i.index("Detection Q Value")
            z = split_i.index("Precursor Charge")
            raw_file = split_i.index("Replicate")
            miss_cleave_idx = split_i.index("Missed Cleavages")
            mz = split_i.index("Precursor Mz")
            rt = split_i.index("Peptide Retention Time")
            return pep_idx, pro_idx, qvalue_idx, z, miss_cleave_idx, mz, rt, raw_file, mod_pep_idx


def ext_skyline_results(infile):
    dicts = {}
    dia_file = {}
    a = get_header_idx(infile)
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            if split_i[a[-1]] + '@' + split_i[a[1]] + '@' + split_i[a[3]] not in dicts:
                dicts[split_i[a[1]] + '@' + split_i[a[-1]] + '@' + split_i[a[3]]] = [split_i]
            else:
                dicts[split_i[a[1]] + '@' + split_i[a[-1]] + '@' + split_i[a[3]]].append(split_i)

            if split_i[a[-2]] not in dia_file:
                dia_file[split_i[a[-2]]] = [split_i]
            else:
                dia_file[split_i[a[-2]]].append(split_i)

    pep = {k.split('@')[1]:k.split('@')[0] + '@' + k.split('@')[2] for k, v in dicts.items() for j in v if j[3] != "#N/A" if float(j[a[2]]) < 0.01}
    pro = {k.split('@')[0]:k.split('@')[0] + '@' + k.split('@')[2] for k, v in dicts.items() for j in v if j[3] != "#N/A" if float(j[a[2]]) < 0.01}
    output = []
    z = {}
    mmc = {}
    output1 = []
    summary = []
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
    length = {}
    for k, v in dia_file.items():
        raw_file.append(k)
        precursors = {}
        for j in v:
            if j[a[0]] + '@' + j[a[-1]] + '@' + j[a[1]] + '@' + j[a[5]] + '@' + j[a[3]] + '@' + j[a[6]] not in precursors:
                precursors[j[a[0]] + '@' + j[a[-1]] + '@' + j[a[1]] + '@' + j[a[5]] + '@' + j[a[3]] + '@' + j[a[6]]] = [j]
            else:
                precursors[j[a[0]] + '@' + j[a[-1]] + '@' + j[a[1]] + '@' + j[a[5]] + '@' + j[a[3]] + '@' + j[a[6]]].append(j)

            if len(j[a[-1]]) not in length:
                length[len(j[a[-1]])] = [j[a[-1]]]
            else:
                length[len(j[a[-1]])].append(j[a[-1]])

        rawfile_pro = {}
        rawfile_pep = {}
        for m, n in precursors.items():
            fdr = {j[a[2]]:j[a[2]] for j in n}
            rawfile_pro[m.split('@')[2]] = 1
            rawfile_pep[m.split('@')[0]] = 1
            output1.append([k] + m.split('@') + [';'.join(list(fdr))])

        print (f"Found {str(len(precursors))} precursors, {len(rawfile_pep)} peptides corresponding to {len(rawfile_pro)} proteins with peptide q-value < 0.01 in the raw file {k}")
        
        summary.append([k, str(len(precursors)), str(len(rawfile_pep)), str(len(rawfile_pro))])
    
    summary.append(['Non-redundant Total', str(len(dicts)), str(len(pep)), str(len(pro))])
    
    write0 = open("{0}_peptide_length_summary.txt".format(infile.rstrip('.tsv')), 'w')
    write0.write("Peptide length" + '\t' + "No. of peptides" + '\n')
    for k, v in length.items():
        #print ("Peptide length " + str(k) + " = " + str(len(v)) + " peptide precursors")
        write0.write(str(k) + '\t' + str(len(v)) + '\n')
    write0.close()  

    write1 = open("{0}_Charge_state_summary.txt".format(infile.rstrip('.tsv')), 'w')
    write1.write("Charge (z)" + '\t' + "Peptide precursors" + '\n')
    for k, v in z.items():
        print ("Charge " + k + " = " + str(len(v)) + " peptide precursors")
        write1.write(k + '\t' + str(len(v)) + '\n')
    write1.close()

    write2 = open("{0}_missed_clavage_summary.txt".format(infile.rstrip('.tsv')), 'w')
    write2.write("Missed cleavage" + '\t' + "Peptide precursors" + '\n')
    for k, v in mmc.items():
        print (f"There are {str(len(v))} peptide precursors with missed cleavage {k}")
        write2.write(k + '\t' + str(len(v)) + '\n')
    write2.close()
                  
    outfile = "{0}_peptides_and_proteins.txt".format(infile.rstrip('.tsv'))
    with open(outfile, 'w') as outf:
        outf.write('Protein\tPeptide\tPrecursor\tCharge\tMissed cleave\tQ-value\n')
        outf.writelines('\t'.join(i) + '\n' for i in output)

    outfile1 = "{0}_raw_file_specific_peptides_and_proteins.txt".format(infile.rstrip('.tsv'))
    with open(outfile1, 'w') as outf1:
        outf1.write('Raw File\tPeptide\tModified Peptide\tProtein\tPrecursor\tCharge(z)\tRT\tQ-values\n')
        outf1.writelines('\t'.join(i) + '\n' for i in output1)

    summary_head = ['Raw File','Peptide Precursors','Peptides','Proteins']
    summrayfile = "{0}_summary.txt".format(infile.rstrip('.tsv'))
    with open(summrayfile, 'w') as sumf:
        sumf.write('\t'.join(summary_head) + '\n')
        sumf.writelines('\t'.join(i) + '\n' for i in summary)

        print (f"In total, there are {summary[-1][1]} peptide precursors {summary[-1][2]} peptide sequences corresponding to {summary[-1][3]} proteins observed")

if __name__== "__main__":
    ext_skyline_results(args.infile[0])
