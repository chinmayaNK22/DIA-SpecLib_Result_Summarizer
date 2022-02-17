from itertools import islice
import argparse

parser = argparse.ArgumentParser(description='''Generate a input matrix for the calculation of protein abundance using iq R scrip''')

parser.add_argument('infile', metavar='-i', type=str, nargs='+', help='Skyline output from DIA-Spectral Library search')
parser.add_argument('conditions', metavar='-c', type=str, nargs='+', help='A .txt file consisting the details of condition specific raw files used for the search')

args = parser.parse_args()

def get_header_idx(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            split_i = i.rstrip().split('\t')
            try:
                mod_pep_idx = split_i.index("Modified Sequence")
                pro_idx = split_i.index("Protein Name")
                transition_idx = split_i.index("Fragment Ion")
                area_idx = split_i.index("Area")
                prec_z = split_i.index("Precursor Charge")
                product_z = split_i.index("Product Charge")
                raw_file = split_i.index("Replicate")
                gene_idx = split_i.index("Protein Gene")
                return raw_file, pro_idx, mod_pep_idx, prec_z, transition_idx, product_z, area_idx, gene_idx
            
            except:
                print ("Some of the essential columns are missing in the exported results from Skyline document")


def get_header(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            print (i.rstrip().split('\t'))

def gen_iq_input(infile, conditions):
    a = get_header_idx(infile)
    dicts = {}
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            sample = split_i[a[0]]
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
                    output.append([j[a[0]], j[a[1]], j[a[2]], j[a[3]], j[a[4]], j[5], j[6], j[7], j[1]])

    print (len(output))

    outfile = "{0}_iq_input.txt".format(infile.rstrip('.txt'))
    with open(outfile, 'w') as outf:
        outf.write('R.Condition\tPG.ProteinGroups\tEG.ModifiedSequence\tFG.Charge\tF.FrgIon\tF.Charge\tF.PeakArea\tPG.Genes\tPG.ProteinNames\n')
        outf.writelines('\t'.join(i) + '\n' for i in output)

if __name__== "__main__":
    gen_iq_input(args.infile[0], args.conditions[0])
