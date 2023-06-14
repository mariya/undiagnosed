import csv
import random

EPIC_DIR="../methylation/files/WE-3687_230602_ResultReport"
EPIC_TSV="%s/WE-3687_230602_SampleMethylationProfile.txt"%EPIC_DIR
OUT_DIR="./out"

def append_bed_line(bed_file, chr, start, end, id, beta_value):
    try:
        if float(start) > 0 and float(end) > 0:
            line = '\t'.join([chr, start, end, id, beta_value])
            with open(bed_file, 'a') as file:
                file.write(line + '\n')
    except:
        return

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
            bed_file = "%s/%s.unsorted.bed" % (OUT_DIR, sample)
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
                    append_bed_line(bed_file, chr, pos, pos, id, beta_value)
                else:
                    # Format: chr18:23503197-23504189
                    # or: chr6:6006456-6006810;chr6:6007154-6007564;chr6:6008624-6009364
                    islands_list = islands.split(';')
                    for isle in islands_list:
                        positions = isle.split(":")[1].split("-")
                        append_bed_line(bed_file, chr, positions[0], positions[1], id, beta_value)

    # Sort all the bedfiles
    print("\n\nSorting BED files...")
    for sample, beta_col in beta_columns.items():
        unsorted_bed_file = "%s/%s.unsorted.bed" % (OUT_DIR, sample)
        sorted_bed_file = "%s/%s.sorted.bed" % (OUT_DIR, sample)

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


