import gzip
import argparse

def convert(filename, outputfilename):
#filename = "./HG002.truvari_collapsed_chr1.vcf.gz"
    if filename.endswith(".gz") or filename.endswith(".bgz"):
        with gzip.open(filename, "r") as fp:
            data = fp.readlines()
    else:
        with open(filename, "r") as fp:
            data = fp.readlines()

    output = []
    for line in data:
        if filename.endswith(".gz") or filename.endswith(".bgz"):
            line = line.decode()
        if line.startswith("#"):
            output.append(line)
            continue
        itemlist = line.split("\t")
        itemlist[3] = itemlist[3].upper()
        itemlist[4] = itemlist[4].upper()
        info = itemlist[7]
        infolist = info.split(';')

        if len(itemlist) < 10:
            continue
        if not itemlist[8].startswith("GT"):
            continue
        # exclude UNK
        for item in infolist:
            if item.startswith("SVTYPE"):
                sv_type = item.split("=")[1]
                break
        else:
            sv_type = None
        if sv_type == "UNK" or sv_type == "CNV":
            continue
        output.append("\t".join(itemlist))

    # outputfilename = "./HG002.truvari_collapsed_converted_chr1.vcf"
    with open(outputfilename, 'w') as fp:
        for line in output:
            fp.write(line)

def parse_args():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Input VCF")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output VCF (stdout)")
    
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    inputfilename = args.input
    outputfilename = args.output
    convert(inputfilename, outputfilename)
    return
if __name__ == '__main__':
    main()
        