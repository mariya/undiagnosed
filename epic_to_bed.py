import csv
import random
from pyliftover import LiftOver

EPIC_DIR="../methylation/files/WE-3687_230602_ResultReport"
EPIC_TSV="%s/WE-3687_230602_SampleMethylationProfile.txt"%EPIC_DIR
OUT_DIR="./out"
LIFTOVER = LiftOver('hg38', 'hg19')

HG19 = "hg19"
HG38 = "hg38"

def bed_filename(sample_name, ref_genome_version, sorted=False):
    sorted = "sorted" if sorted else "unsorted"
    return "%s/%s.%s.%s.bed" % (OUT_DIR, sample, ref_genome_version, sorted)
    
def to_hg19(chr, pos):
    converted = LIFTOVER.convert_coordinate(chr, pos)
    try:
        return converted[0][1]
    except:
        print("Could not find hg19 coordinates for hg38 coordinates %s:%s" % (chr, pos))
        return 0

def append_bed_line(sample, chr, start, end, id, beta_value):
    bed_file_hg38 = bed_filename(sample, HG38)
    bed_file_hg19 = bed_filename(sample, HG19)
    
    try:
        i_start = int(start)
        i_end = int(end)
        if i_start > 0 and i_end > 0:
            # hg38
            line = '\t'.join([chr, start, end, id, beta_value])
            with open(bed_file_hg38, 'a') as file:
                file.write(line + '\n')

            # hg19 conversion
            hg19_start = to_hg19(chr, i_start)
            hg19_end = to_hg19(chr, i_end)
            if hg19_start > 0 and hg19_end > 0:
                line = '\t'.join([chr, f"{hg19_start}", f"{hg19_end}", id, beta_value])
                with open(bed_file_hg19, 'a') as file:
                    file.write(line + '\n')
    except ValueError:
        return

def sort_bed_file(sample, ref_genome_version):
    unsorted_bed_file = bed_filename(sample, ref_genome_version, False)
    sorted_bed_file = bed_filename(sample, ref_genome_version, True)

    # Add an IGV track line to each bed file, in a random color
    r = "%s" % random.randint(0,255)
    g = "%s" % random.randint(0,255)
    b = "%s" % random.randint(0,255)
    rgb = ','.join([r,g,b])

    track_line = 'track name="%s EPIC" description="%s EPIC" graphType="bar" color=%s viewLimits=0:1 maxHeightPixels=50:50:50' % (sample, sample, rgb)
    with open(sorted_bed_file, 'w') as file:
        file.write(track_line + '\n')

    with open(unsorted_bed_file, 'r') as file:
        lines = file.readlines()
        sorted_lines = sorted(lines, key=lambda x: (x.split('\t')[0], int(x.split('\t')[1])))
        with open(sorted_bed_file, 'a') as file:
            file.writelines(sorted_lines)


with open(EPIC_TSV, newline='') as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter='\t')

    # Get header row
    header = next(tsvreader)

    # Get the relevant sample names
    print("Extracting column names...")
    beta_columns = {}
    for index, value in enumerate(header):
        str="%i %s" % (index, value)
        print(str)
        if "AVG_Beta" in value:
            sample_name = value.split(".")[0]
            beta_columns[sample_name] = index

    print("\n\nGenerating BED files per sample...")
    # Loop through all rows
    # Columns of interest:
    # - TargetID (1)
    # - [SAMPLE].AVG_Beta - lookup in beta_columns
    # - CHR (163)
    # - MAPINFO (164)
    # - UCSC_CPG_ISLANDS_NAME (174)

    for row in tsvreader:
        id = row[1]
        chr = row[163]
        pos = row[164]
        islands = row[174]

        # Loop through all samples
        for sample, beta_col in beta_columns.items():
            beta_value = row[beta_col]

            # Discard any rows with invalid beta value
            try:
                beta_float = float(beta_value)
            except:
                beta_float = -1

            if beta_float >= 0 and beta_float <= 1:
                # UCSC_CPG_ISLANDS_NAME
                # If this field is blank, use pos for start and end
                if islands.strip() == "":
                    append_bed_line(sample, chr, pos, pos, id, beta_value)
                else:
                    # Format: chr18:23503197-23504189
                    # or: chr6:6006456-6006810;chr6:6007154-6007564;chr6:6008624-6009364
                    islands_list = islands.split(';')
                    for isle in islands_list:
                        positions = isle.split(":")[1].split("-")
                        append_bed_line(sample, chr, positions[0], positions[1], id, beta_value)

    # Sort all the bedfiles
    print("\n\nSorting BED files...")
    for sample, beta_col in beta_columns.items():
        sort_bed_file(sample, HG38)
        sort_bed_file(sample, HG19)

