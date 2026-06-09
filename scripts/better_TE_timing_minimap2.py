from collections import defaultdict
from pathlib import Path
import re
import sys
import os
import math
import tempfile
import subprocess
import time
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from shutil import copyfile

def argsParser():
    import argparse
    parser = argparse.ArgumentParser(description='Calculates Jukes-Cantor and Gap-compressed identity\ngiven gff of TEs in EarlGrey format')
    parser.add_argument('-i', '--in_gff', type=str, required=True,
                        help='Input gff, nested TEs fixed')
    parser.add_argument('-o', '--out_gff', type=str, required=True,
                        help='Output gff)')
    parser.add_argument('-g', '--genome', type=str, required=True,
                        help='path to genome assembly)')
    parser.add_argument('-d', '--debug', action="store_true",
                        help='Save temporary files for debugging to tmp')
    parser.add_argument('-td', '--tmp_dir', type=str, default="tmp/",
                        help='If debugging, path directory for temporary files)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='number of threads to run minimap2 with')
    parser.add_argument('--k', '-k', type=int, default=10,
                        help='k value for minimap2')
    parser.add_argument('--w', '-w', type=int, default=5,
                        help='w value for minimap2')
    parser.add_argument('--m', '-m', type=int, default=5,
                        help='m value for minimap2')
    parser.add_argument('-c', '--coverage', type=float, default=0.8,
                        help='Coverage queries must have by target to calculate divergence')
    args = parser.parse_args()

    if args.coverage <= 0 or args.coverage > 1:
        sys.exit("Coverage must be between 0 and 1")
    if args.in_gff == args.out_gff:
        sys.exit("In and out gffs must be different paths")

    return(args)

def make_fastas(gff_file, genome_path, tmp_dir, assessable):
    # FASTA index for efficient random access
    print(f"Creating fasta")
    genome = Fasta(genome_path, sequence_always_upper=True)

    # Store sequences by repeat type
    records_by_type = defaultdict(list)

    with open(gff_file) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            
            # get coordindates, type etc from gff
            chrom, source, repeat_type, start, end, score, strand, phase, attrs = fields[:9]
            start = int(start)
            end = int(end)
            repeat_subclass = re.sub('/.*', '', repeat_type)

            # Only assess TEs (ignore simple etc)
            if repeat_subclass not in assessable:
                continue
            
            # To avoid including tandam arrays, check if repeat length <= 2* consensus length
            metadata = attrs.split(";")
            con_len = [s for s in metadata if "CON_LEN" in s]
            con_len = int(re.sub("CON_LEN=", "", *con_len))
            if end - start >= con_len :
                continue

            # GFF uses 1-based inclusive coordinates
            seq = genome[chrom][start - 1:end].seq

            # Reverse-complement if feature is on minus strand
            if strand == "-":
                seq = str(Seq(seq).reverse_complement())

            # Create a FASTA record
            rec_id = f"{chrom}:{start}-{end}:{strand}"
            record = SeqRecord(Seq(seq), id=rec_id, description="")

            # Group by the GFF type column
            records_by_type[repeat_type].append(record)
    repeat_types = []

    # Write one multi-FASTA file per repeat type
    for repeat_type, records in records_by_type.items():
        safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", repeat_type) # Replace backslash in repeat type with underscore
        out_fasta = f"{tmp_dir}/{safe_name}.fa"
        SeqIO.write(records, out_fasta, "fasta")
        repeat_types.append(safe_name)
    return(repeat_types)

def calc_matches_mismatches(nm, cigar_string):
    # Regex to find all occurrences of a number followed by 'I' or 'D', or a number folowed by an 'M'
    # Count matches
    match = r'(\d+)[M]'
    matches = re.findall(match, cigar_string, re.IGNORECASE)
    matches_bp = sum(int(num) for num in matches)
    # Count insertions and deletions from cigar string
    if 'I' in cigar_string or 'D' in cigar_string:
        indel = r'(\d+)[ID]'
        indels = re.findall(indel, cigar_string, re.IGNORECASE)
        indels_bp = sum(int(num) for num in indels)
    else:
        indels_bp = 0
    
    # NM = I + D + Mismatches
    mismatches_bp = nm - indels_bp
    return matches_bp, mismatches_bp

def jukes_cantor(matches, mismatches):
    p = mismatches / (matches + mismatches)
    return -0.75 * math.log(1.0 - (4.0 / 3.0) * p)

def align_calc(repeat_types, tmp_dir, threads, coverage, k, w, m):
    # create dictionary for best hits
    best_hits = {}

    # loop through repeat familys, aligning with minimap2 and then calculating best hit
    for repeat_type in repeat_types:
        print(f"Running calculations for {repeat_type}")
        # get path to superfamily fasta
        fasta_path=f"{tmp_dir}/{repeat_type}.fa"
        paf_path=f"{tmp_dir}/{repeat_type}.paf"

        # Runs all-vs-all minimap2 alignment for repeats in superfamily
        cmd = [
            "minimap2",
            "-c", "-P", "-D", "-k", str(k), "-w", str(w), "-m", str(m), "-t", str(threads),
            fasta_path,
            fasta_path
        ]
        
        with open(paf_path, "w") as outpaf:
                subprocess.run(cmd, stdout=outpaf, stderr = subprocess.DEVNULL, check=True)

        # Read through paf and calculate jc distance
        with open(paf_path) as fh:
                for line in fh:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 12:
                        continue
                    # pull out query name, start, end and length
                    qname=cols[0]
                    qlen=int(cols[1])
                    qstart=int(cols[2])
                    qend=int(cols[3])
                    tname=cols[5]
                    
                    # If already calculate skip
                    if qname in best_hits.keys():
                        continue

                    # skip to next alignment if < specified coverage (default 80%)
                    if qend - qstart < coverage * qlen:
                        continue

                    # Parse extra tags and pull out cigar string and nm
                    tags = {}
                    for field in cols[12:]:
                        parts = field.split(":", 2)
                        if len(parts) == 3:
                            tags[parts[0]] = parts[2]

                    # pull out minimap2 cigar string and NM to calc jc_dist
                    cg = tags.get("cg")
                    nm = int(tags.get("NM"))
                    matches, mismatches = calc_matches_mismatches(nm, cg)
                    jc_dist = jukes_cantor(matches, mismatches)
                    # pull out minimap2 calculated distance ("Gap-compressed identity")
                    de = float(tags.get("de"))
                    
                    best_hits[qname] = [str(abs(jc_dist)), str(de)]
    return(best_hits)

def write_gff(gff_file, out_gff_path, best_hits):
    print(f"Writing gff")
    with open(gff_file, "r") as fh, open(out_gff_path, "w") as o:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            line = line.rstrip("\n")

            fields = line.split("\t")
            if len(fields) < 9:
                continue
            
            # get coordindates, type etc from gff
            chrom, source, repeat_type, start, end, score, strand, phase, attrs = fields[:9]
            start = int(start)
            end = int(end)
            fields[1] = "EarlGrey"
            
            repeat=f"{chrom}:{start}-{end}:{strand}"
            
            # get JC distance and Gap-compressed identity from best hits if present
            joined_fields="\t".join(fields)
            if repeat in best_hits:
                jc_dist, de = float(best_hits[repeat][0]), float(best_hits[repeat][1])
                out_line=f"{joined_fields};JC_DIST={jc_dist:.3f};GC_ID={de:.3f}\n"
            else:
                out_line=f"{joined_fields}\n"
            o.write(out_line)

def main():
    
    args=argsParser()
    
    # Make temporary directory, if debugging use desired path
    if args.debug:
        tmp_dir=args.tmp_dir
        if not os.path.exists(args.tmp_dir):
            os.makedirs(args.tmp_dir)
    else:
        tmp = tempfile.TemporaryDirectory(delete=False)
        tmp_dir = tmp.name
    if args.debug:
        print(tmp_dir)
    # Set repeat subclasses to be assesed
    assessable = ['LINE', 'LTR', 'RC', 'DNA', 'SINE', 'PLE', 'Unknown']
    # make fasta for each repeat type
    repeat_types = make_fastas(args.in_gff, args.genome, tmp_dir, assessable)
    # align repeats to each other and calculate JC_dist
    start_time = time.time()
    best_hits = align_calc(repeat_types, tmp_dir, args.threads, args.coverage, str(args.k), str(args.w), str(args.m))
    end_time = time.time()
    if args.debug:
        print(f"Alignment took {end_time - start_time}seconds")

    # Add newly calculated JC_dist to gff
    write_gff(args.in_gff, args.out_gff, best_hits)

    # cleanup temporary directory if not debugging
    if not args.debug:
        tmp.cleanup()

if __name__ == '__main__':
    __version__ = '0.1'
    try:
        main()
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user\n")
        exit(-1)