import argparse
import json
import re

parser = argparse.ArgumentParser(description='Compose a samplesheet given sample metadata.')
parser.add_argument('-i', '--input', help="Input json file with sample metadata", default="0")
parser.add_argument('-m', '--mask', help="Base Mask", required=False)
parser.add_argument('-l', '--lanes', help="Lanes", required=False)
args = parser.parse_args()


'''
   Handle both types of barcodes, dual or single
'''
lanes = []
haveLanes = ""

with open("samplesheet.csv", "w") as ss:
    ss.write("[Data]\n")
    ss_lines = []
    lanes = re.split(",", args.lanes) if args.lanes else []
    ss_fields = 0
    dualBarcodes = False

    with open(args.input) as js:
        sample = json.load(js)
        print(sample)
        name = sample['name']
        barcodes = sample['barcodes']
        for barcode in barcodes:
            try:
                (bc1, bc2) = barcode.split("-")
                ss_lines.append(f'{name},{bc1},{bc2}\n')
                dualBarcodes = True
            except:
                ss_lines.append(f'{name},{barcode}\n')

    if lanes and len(lanes) > 0:
        haveLanes = "Lane,"
        ss_fields += 1
    if dualBarcodes:
        ss.write(haveLanes + "Sample_ID,index,index2\n")
        ss_fields += 3
    else:
        ss.write(haveLanes + "Sample_ID,index\n")
        ss_fields += 2

    if lanes and len(lanes) > 0:
        for lane in lanes:
            for line in ss_lines:
                out_line = ",".join([lane, line])
                ss.write(out_line)
    else:
        for line in ss_lines:
            ss.write(line)

    if args.mask and re.match('Y', args.mask):
        print(f'We have a valid mask {args.mask}')
        mask_line = ",".join(["OverrideCycles", args.mask])
        for i in range(0, ss_fields - 2):
            mask_line += ","
        ss.write("\n\n[Settings]\n")
        ss.write(mask_line + "\n")
    ss.close()
