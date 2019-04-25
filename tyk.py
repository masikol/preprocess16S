import os

# outdir_path = "/home/deynonih/Documents/Univier/Courseache/preprocess16S_result_2019-04-11_11.27.55"
outdir_path = "/home/deynonih/Documents/Univier/Courseache/Metagenomics/MiSeq_16S_Test_data/FASTQ_Generation/ATCCMiSeq2x300R10_L001-ds.63c79c7933b147f784741b8964048871"

names = dict()
# names["R1"] = "08_MG-fecal-A2_S8_L001_R1_001"
# names["R2"] = "08_MG-fecal-A2_S8_L001_R2_001"
names["R1"] = "ATCCMiSeq2x300R10_S10_L001_R1_001"
names["R2"] = "ATCCMiSeq2x300R10_S10_L001_R2_001"
# names["R1"] = "fecal_trimmed_R1_paired"
# names["R2"] = "fecal_trimmed_R2_paired"


# result_paths = {
#     "mR1": "{}{}{}.16S.fastq".format(outdir_path, os.sep, names["R1"]),
#     "trR1": "{}{}{}.trash.fastq".format(outdir_path, os.sep, names["R1"]),
#     "mR2": "{}{}{}.16S.fastq".format(outdir_path, os.sep, names["R2"]),
#     "trR2": "{}{}{}.trash.fastq".format(outdir_path, os.sep, names["R2"]),
# }

result_paths = {
    "mR1": "{}{}{}.fastq".format(outdir_path, os.sep, names["R1"]),
    "trR1": "{}{}{}.trash.fastq".format(outdir_path, os.sep, names["R1"]),
    "mR2": "{}{}{}.fastq".format(outdir_path, os.sep, names["R2"]),
    "trR2": "{}{}{}.trash.fastq".format(outdir_path, os.sep, names["R2"]),
}

# num_reads = int(456802 / 2)
num_reads = 777620
#num_reads = 31403



# --------------------------------------------------------------------------------------

filt_read_paths = {
    "mR1": result_paths["mR1"],
    "mR2": result_paths["mR2"]
}
filt_read_files = dict()
try:
    filt_read_files["mR1"] = open(filt_read_paths["mR1"], 'r')
    filt_read_files["mR2"] = open(filt_read_paths["mR2"], 'r')
except OSError as oserror:
    print("Error while opening one of result files", repr(oserror))
    print("Here ary your files:\n\t{}\n\t{}".format(filt_read_paths["mR1"], filt_read_paths["mR2"]))
    for key in filt_read_files.keys():
        filt_read_files[key].close()
    input("Press enter to exit:")
    exit(1)

more_common_name = names["R1"][: names["R1"].find("_R1_")]
merged_path = "{}{}{}.merged.fastq".format(outdir_path, os.sep, more_common_name)
blast_rep = "blastn_report.txt"
fasta_rep = "fasta36_report.txt"
query = "query.fasta"
sbjct = "subject.fasta"
ncbi_fmt_db = "SILVA_db/SILVA_132_SSURef_Nr99_tax_silva.fasta"
constV3V4_path = "/home/deynonih/Documents/Univier/Courseache/preprocess16S_result_2019-04-11_11.27.55/constant_region_V3-V4.fasta"

cmd_for_blastn = """blastn -query {} -db {} -penalty -1 -reward 2 -ungapped \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sacc sstrand" \
-out {} -task blastn -max_target_seqs 1""".format(query, ncbi_fmt_db, blast_rep)
QSEQID, SSEQID, PIDENT, LENGTH, MISMATCH, GAPOPEN, QSTART, QEND, SSTART, SEND, EVALUE, BITSCORE, SACC, SSTRAND = range(14)
draft_cmd_for_blastdbcmd = "blastdbcmd -db SILVA_db/SILVA_132_SSURef_Nr99_tax_silva.fasta -entry REPLACE_ME -out {}".format(sbjct)
cmd_for_fasta = "fasta36 {} {} -n -r +4/-2 -f 20 -g 10 -m 8 -3 > {}".format(query, sbjct, fasta_rep)
# cmd_for_fasta = "fasta36 {} {} -n -f 20 -g 10 -3 -O {}".format(query, sbjct, fasta_rep)


RC_DICT = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'U': 'A',   # in this case it does not matter
    'N': 'N'    # just in case
}
single_nucl_rc = lambda nucl: RC_DICT[nucl]
rc = lambda seq: "".join(map(single_nucl_rc, seq[::-1]))


const_V3_V4 = ""
with open(constV3V4_path, 'r') as const_seq_file:
    const_seq_file.readline()
    for line in const_seq_file:
        const_V3_V4 += line.strip()
constV3V4_len = len(const_V3_V4)
del const_V3_V4


MAX_ALIGN_OFFSET = 40
MIN_OVERLAP = 10
MAX_CONST_MID_OFFS = 70
MIN_CONST_COV = 0.90
MIN_CONST_IDENT = 80.0



# ==================================================================================================
# functions


def merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual):
    # print('\n' + '=' * 60 + '\n')
    print(fseq[loffset : loffset + overl] + '\n')
    print(rseq[:overl])
    print('\n' + '=' * 60)

    # beginning comes from forward read
    merged_seq = fseq[: loffset]
    merged_qual = fqual[: loffset]

    # "recover" overlapping region depending on quality 
    for i, j in zip(range(loffset, loffset + overl), range(overl)):
        try:
            qual, nucl = (fqual[i], fseq[i]) if ord(fqual[i]) >= ord(rqual[j]) else (rqual[j], rseq[j])
        except IndexError:
            print(i, j)
            print(loffset, overl)
            print(len(fseq))

            with open(query, 'w') as query_file:
                query_file.write(fastq_recs["mR1"]["seq_id"].replace('@', '>') + '\n')
                query_file.write(fseq)
            with open(sbjct, 'w') as sbjct_file:
                sbjct_file.write(fastq_recs["mR2"]["seq_id"].replace('@', '>') + '\n')
                sbjct_file.write(rseq)
            os.system("fasta36 {} {} -n -r +5/-2 -f 20 -g 10 -3 > {}".format(query, sbjct, fasta_rep))

            exit(0)
        merged_seq += nucl
        merged_qual += qual

    # end comes from reverse read
    merged_seq += rseq[overl :]
    merged_qual += rqual[overl :]

    return (merged_seq, merged_qual)


def blast_and_align(rseq):
    """
    Assume that query is still the same
    """

    # blast forward read 
    os.system(cmd_for_blastn)   # query is still the same
    with open(blast_rep, 'r') as blast_repf:
        faref_report = blast_repf.readline().strip().split('\t')
        # print('\n' + '~' * 40 + '\n')
        # print(faref_report)

        # get ref sequence
        cmd_for_blastbdcmd = draft_cmd_for_blastdbcmd.replace("REPLACE_ME", faref_report[SACC])
        os.system(cmd_for_blastbdcmd)
        if faref_report[SSTRAND] == "minus":
            with open(sbjct, 'r') as sbjct_file:
                sbjct_seq_id = sbjct_file.readline().strip()
                sbjct_seq = ""
                for line in sbjct_file:
                    sbjct_seq += line.strip()
            sbjct_seq = rc(sbjct_seq)
            with open(sbjct, 'w') as sbjct_file:
                sbjct_file.write(sbjct_seq_id + '\n')
                sbjct_file.write(sbjct_seq)

        # align reverse read against the reference
        with open(query, 'w') as query_file:
            query_file.write(fastq_recs["mR2"]["seq_id"].replace('@', '>') + '\n')
            query_file.write(rseq)
        os.system(cmd_for_fasta)

        with open(fasta_rep, 'r') as fasta_repf:
             # raref means reverse [read] against reference [sequence]
            raref_report = fasta_repf.readline().strip().split('\t')

        return faref_report, raref_report


def search_for_constant_region(loffset, overl, fseq, fqual, rseq, rqual):

    # try to merge
    merged_seq, merged_qual = merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual)

    # form query file
    with open(sbjct, 'w') as sbjct_file:
        sbjct_file.write(">presumably_merged_read\n")
        sbjct_file.write(merged_seq)

    # form subject file
    os.system("cat {} > {}".format(constV3V4_path, query))

    # align
    os.system(cmd_for_fasta)

    # parse report
    with open(fasta_rep, 'r') as fasta_repf:
        pmer_report = fasta_repf.readline().strip().split('\t') # "pmer" means Presumably MERged [read]

    # check if there are a constant region
    enough_cov = float(pmer_report[LENGTH]) / constV3V4_len > MIN_CONST_COV
    enough_indent = float(pmer_report[PIDENT]) > MIN_CONST_IDENT
    in_the_midl = int(pmer_report[SSTART]) > MAX_CONST_MID_OFFS
    in_the_midr = int(pmer_report[SEND]) > len(merged_seq) - MAX_CONST_MID_OFFS

    if enough_cov and enough_indent and in_the_midl and in_the_midr:
        print("\nConstant region found\n")
        # WRITE FASTQ RECORD
    else:
        print("\nPUTATIVE CHIMERA!!\n")
        print("R1:  " + fseq)
        print("R2:  " + rseq)

        # WRITE FASTQ RECORD



# process

with open(merged_path, 'w') as merged_file:
    # for key in filt_read_files.keys():
    #     for i in range(4 * ):
    #         filt_read_files[key].readline()
    reads_processed = 0
    while reads_processed < num_reads:
        print(reads_processed)
        fastq_recs = dict()           # this dict should consist of two fastq-records: from R1 and from R2
        for key in filt_read_files.keys():
            fastq_recs[key] = {                    #read all 4 lines of fastq-record
                "seq_id": filt_read_files[key].readline().strip(),
                "seq": filt_read_files[key].readline().strip(),
                "optional_id": filt_read_files[key].readline().strip(),
                "quality_str": filt_read_files[key].readline().strip()
            }

        fseq = fastq_recs["mR1"]["seq"]
        fqual = fastq_recs["mR1"]["quality_str"]
        rseq = rc(fastq_recs["mR2"]["seq"])
        rqual = fastq_recs["mR2"]["quality_str"][::-1]

        # Firstly align forward read against reverse.

        # form query file
        with open(query, 'w') as query_file:
            query_file.write(fastq_recs["mR1"]["seq_id"].replace('@', '>') + '\n')
            query_file.write(fseq)

        # form subject file
        with open(sbjct, 'w') as sbjct_file:
            sbjct_file.write(fastq_recs["mR2"]["seq_id"].replace('@', '>') + '\n')
            sbjct_file.write(rseq)

        # align forward read against reverse read
        os.system(cmd_for_fasta)

        # parse report
        with open(fasta_rep, 'r') as fasta_repf:
            far_report = fasta_repf.readline().strip().split('\t')  # "far" means Forward [read] Against Reverse [read]


        # Check how they have aligned
        # naive assumption
        reads_normally_overlap = True   
        # discard following situation:
        # --FFFFFF
        # RRRRRR--
        too_short = int(far_report[QSTART]) < int(far_report[SSTART])
        too_short = too_short or int(far_report[QEND]) / len(fseq) < int(far_report[SEND]) / len(rseq)
        reads_normally_overlap = reads_normally_overlap and not too_short
        # catch randomly occured alignment in center of sequences, 
        # e. i. reads don't align against one another in a proper way
        rand_align = len(fseq) - int(far_report[QEND]) > MAX_ALIGN_OFFSET or int(far_report[SSTART]) > MAX_ALIGN_OFFSET
        # resulting conslusion
        reads_normally_overlap = reads_normally_overlap and not rand_align


        if reads_normally_overlap:
            """(ok, reads overlap where they should do it)
            FFFFFFFF----
            ----RRRRRRRR
            """
            print("\nOVERLAP\n")
            loffset = int(far_report[QSTART]) - int(far_report[SSTART])
            # in the next line we have 1-based terms, hense I need to substract 1
            overl = int(far_report[LENGTH]) + int(far_report[SSTART]) + (len(fseq) - int(far_report[QEND])) - 1
            merged_seq, merged_qual = merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual)
            print('\n' + merged_seq + '\n')

        elif rand_align:
            # (randomly occured alignment in center of sequences, need to blast)
            # "faref" means Forward [read] Against REFerence [sequence]
            # "raref" means Reverse [read] Against REFerence [sequence]
            gap = False
            faref_report, raref_report = blast_and_align(rseq)
            forw_start = int(faref_report[SSTART]) - int(faref_report[QSTART])
            forw_end = forw_start + len(fseq)
            rev_start = int(raref_report[SSTART]) - int(raref_report[QSTART])

            if forw_end < rev_start:
                gap = True

            if not gap and forw_end - rev_start > MIN_OVERLAP:
                # search for constant region
                if forw_start < rev_start:
                    loffset = rev_start - forw_start
                    overl = forw_end - rev_start 
                    search_for_constant_region(loffset, overl, fseq, fqual, rseq, rqual)
                else:
                    print("\nPUTATIVE CHIMERA!!!\n")
                    print("R1:  " + fseq)
                    print("R2:  " + rseq)

                    # WRITE FASTQ RECORD

            elif not gap and forw_end - rev_start <= MIN_OVERLAP:
                # length of the overlapping region is small,
                # but reverse read alignes as they should do it,
                # so let it be
                print("Let it be\n")
                loffset = rev_start - forw_start
                overl = forw_end - rev_start 
                merged_seq, merged_qual = merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual)
                # WRITE FASTQ RECORD
                
            elif gap:
                # here we have a gap
                print("\nGAP\n")
                # fill in the gap with N
                gap_len = rev_start - forw_end
                merged_seq = fseq + 'N' * (gap_len - 1) + rseq
                # Illumina uses Phred33
                merged_qual = fqual + chr(33) * (gap_len - 1) + rqual
                print(faref_report)
                print(raref_report)
                print(merged_seq)
            else:
                print("Unforeseen case! Report it to me -- I will fix it.")
                exit(0)

        elif too_short:
            print("\nTOO SHORT!")
            """(too short)
            --FFFFFFF
            RRRRRRR--
            """
            # WRITE FASTQ RECORD
        else:
            print("Unforeseen case! Report it to me -- I will fix it.")
            exit(0)

        if not reads_normally_overlap:
            with open(query, 'w') as query_file:
                query_file.write(fastq_recs["mR1"]["seq_id"].replace('@', '>') + '\n')
                query_file.write(fseq)
            with open(sbjct, 'w') as sbjct_file:
                sbjct_file.write(fastq_recs["mR2"]["seq_id"].replace('@', '>') + '\n')
                sbjct_file.write(rseq)
            os.system("fasta36 {} {} -r +4/-2 -n -f 20 -g 10 -3 > {}".format(query, sbjct, fasta_rep))
            input("Press ENTER:")

        reads_processed += 1
        # input("Press ENTER")
        # print('*' * 60 + '\n\n')



    # merged_file.write(fastq_recs["mR1"]["seq_id"] + '\n')
    # merged_file.write(merged_seq + '\n')
    # merged_file.write(fastq_recs["mR1"]["optional_id"] + '\n')
    # merged_file.write(merged_qual + '\n')


for key in filt_read_files.keys():
    filt_read_files[key].close()