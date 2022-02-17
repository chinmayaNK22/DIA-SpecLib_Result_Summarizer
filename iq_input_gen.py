from itertools import islice

infile = "Reanalysi_02122021\\NMO_DIA_25Da_PI_6_3_0.02Da_20ppm_spec_lib_search_120221_Result.tsv"

conditions = "NMO_Sample-Condition_Info.txt"

def get_header(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            print (i.rstrip().split('\t'))

get_header(infile)
dicts = {}
with open(infile) as file:
    for i in islice(file, 1, None):
        split_i = i.rstrip().split('\t')
        sample = split_i[-1].split('_')[0]
        if sample not in dicts:
            dicts[sample] = [split_i]
        else:
            dicts[sample].append(split_i)


output = []
with open(conditions) as file:
    for i in islice(file, 1, None):
        split_i = i.rstrip().split('\t')
        if split_i[0] in dicts:
            for j in dicts[split_i[0]]:
                #print (split_i[0], j[0], j[9], j[6], j[13], j[14], j[-5], j[2], j[0])
                output.append([j[-1], j[0], j[9], j[6], j[13], j[14], j[-5], j[2], j[0]])

print (len(output))

outfile = "{0}_iq_input.txt".format(infile.rstrip('.txt'))
with open(outfile, 'w') as outf:
    outf.write('R.Condition\tPG.ProteinGroups\tEG.ModifiedSequence\tFG.Charge\tF.FrgIon\tF.Charge\tF.PeakArea\tPG.Genes\tPG.ProteinNames\n')
    outf.writelines('\t'.join(i) + '\n' for i in output)
