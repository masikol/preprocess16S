import os

outdir_path = "/home/deynonih/Documents/Univier/Courseache/preprocess16S_result_2019-04-11_11.27.55"
# outdir_path = "/home/deynonih/Documents/Univier/Courseache/Metagenomics/MiSeq_16S_Test_data/FASTQ_Generation/ATCCMiSeq2x300R10_L001-ds.63c79c7933b147f784741b8964048871"

names = dict()
names["R1"] = "08_MG-fecal-A2_S8_L001_R1_001"
names["R2"] = "08_MG-fecal-A2_S8_L001_R2_001"
# names["R1"] = "ATCCMiSeq2x300R10_S10_L001_R1_001"
# names["R2"] = "ATCCMiSeq2x300R10_S10_L001_R2_001"

result_paths = {
    "mR1": "{}{}{}.16S.fastq".format(outdir_path, os.sep, names["R1"]),
    "trR1": "{}{}{}.trash.fastq".format(outdir_path, os.sep, names["R1"]),
    "mR2": "{}{}{}.16S.fastq".format(outdir_path, os.sep, names["R2"]),
    "trR2": "{}{}{}.trash.fastq".format(outdir_path, os.sep, names["R2"]),
}

# result_paths = {
#     "mR1": "{}{}{}.fastq".format(outdir_path, os.sep, names["R1"]),
#     "trR1": "{}{}{}.trash.fastq".format(outdir_path, os.sep, names["R1"]),
#     "mR2": "{}{}{}.fastq".format(outdir_path, os.sep, names["R2"]),
#     "trR2": "{}{}{}.trash.fastq".format(outdir_path, os.sep, names["R2"]),
# }

num_reads = int(456802 / 2)

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
fasta_fmt_db = "SILVA_db/SILVA_132_SSURef_Nr99_tax_silva.fasta"

cmd_for_blastn = """blastn -query {} -db {} -penalty -1 -reward 2 -ungapped \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sacc sstrand" \
-out {} -task blastn -max_target_seqs 1""".format(query, ncbi_fmt_db, blast_rep)
QSEQID, SSEQID, PIDENT, LENGTH, MISMATCH, GAPOPEN, QSTART, QEND, SSTART, SEND, EVALUE, BITSCORE, SACC, SSTRAND = range(14)
draft_cmd_for_blastdbcmd = "blastdbcmd -db SILVA_db/SILVA_132_SSURef_Nr99_tax_silva.fasta -entry REPLACE_ME -out {}".format(sbjct)
cmd_for_fasta = "fasta36 {} {} -n -f 20 -g 10 -m 8 -3 -O {}".format(query, sbjct, fasta_rep)
# cmd_for_fasta = "fasta36 {} {} -n -f 20 -g 10 -3 -O {}".format(query, sbjct, fasta_rep)


RC_DICT = {
	'A': 'T',
	'T': 'A',
	'C': 'G',
	'G': 'C',
	'U': 'A',
	'N': 'N'
}
single_nucl_rc = lambda nucl: RC_DICT[nucl]
rc = lambda seq: "".join(map(single_nucl_rc, seq[::-1]))


def ok_overlap(loffset, overl, fseq, fqual, rseq, rqual):
	print('\n' + '=' * 60 + '\n')
	print(fseq[loffset : loffset + overl] + '\n')
	print(rseq[:overl])
	print('\n' + '=' * 60 + '\n')

	merged_seq = fseq[: loffset]
	merged_qual = fqual[: loffset]
	for i, j in zip(range(loffset, loffset + overl), range(overl)):
		qual, nucl = (fqual[i], fseq[i]) if ord(fqual[i]) >= ord(rqual[j]) else (rqual[j], rseq[j])
		merged_seq += nucl
		merged_qual += qual
	merged_seq += rseq[overl :]
	merged_qual += rqual[overl :]
	return (merged_seq, merged_qual)


def align_reads_ag_some_ref(rseq):
	"""
	Assume that query is still the same
	"""

	# blast forward read 
	os.system(cmd_for_blastn)	# query is the same
	with open(blast_rep, 'r') as blast_repf:
		for line in blast_repf:
			# faref_report = blast_repf.readline().strip().split('\t')   # faref means forward [read] against reference [sequence]
			faref_report = line.strip().split('\t')
			print('\n' + '~' * 40 + '\n')
			print(faref_report)

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
				raref_report = fasta_repf.readline().strip().split('\t')	# raref means reverse [read] against reference [sequence]

			print(raref_report)

			print(int(raref_report[LENGTH]) / len(rseq))

			if int(raref_report[LENGTH]) / len(rseq) > 0.7:
				yield faref_report, raref_report

		yield "NOTHING", "MORE"


with open(merged_path, 'w') as merged_file:
	try:
		for key in filt_read_files.keys():
			for i in range(4 * 40):
				filt_read_files[key].readline()


		fastq_recs = dict()           # this dict should consist of two fastq-records: from R1 and from R2
		for key in filt_read_files.keys():
			fastq_recs[key] = { 	               #read all 4 lines of fastq-record
                "seq_id": filt_read_files[key].readline().strip(),
                "seq": filt_read_files[key].readline().strip(),
                "optional_id": filt_read_files[key].readline().strip(),
                "quality_str": filt_read_files[key].readline().strip()
            }
	except IOError as ioerror:
		print("Error while parsing one of the .fastq read files", repr(ioerror))
		for key in filt_read_files.keys():
			filt_read_files[key].close()
		input("Press enter to exit:")
		exit(1)

	fseq = fastq_recs["mR1"]["seq"]
	fqual = fastq_recs["mR1"]["quality_str"]
	rseq = rc(fastq_recs["mR2"]["seq"])
	rqual = fastq_recs["mR2"]["quality_str"][::-1]

	with open(query, 'w') as query_file:
		query_file.write(fastq_recs["mR1"]["seq_id"].replace('@', '>') + '\n')
		query_file.write(fseq)

	with open(sbjct, 'w') as sbjct_file:
		sbjct_file.write(fastq_recs["mR2"]["seq_id"].replace('@', '>') + '\n')
		sbjct_file.write(rseq)

	os.system(cmd_for_fasta)
	with open(fasta_rep, 'r') as fasta_repf:
		far_report = fasta_repf.readline().strip().split('\t')  # far means forward [read] against reverse [read]

	print(far_report)
	max_align_offset = 10

	reads_normally_overlap = True 	# naive assumption
	reads_normally_overlap = reads_normally_overlap and int(far_report[QSTART]) > int(far_report[SSTART])
	complicated_term = len(fseq) - int(far_report[QEND]) > max_align_offset or int(far_report[SSTART]) > max_align_offset
	reads_normally_overlap = reads_normally_overlap and not complicated_term

	if reads_normally_overlap:
		"""(ok, reads overlap)
		FFFFFFFF
			RRRRRRRR
		"""
		print("\nOVERLAP")
		loffset = int(far_report[QSTART]) - int(far_report[SSTART])
		# in the next line we have four 1-based terms, hense I need to substract 4
		overl = int(far_report[SEND]) + int(far_report[SSTART]) + (len(fseq) - int(far_report[QEND])) - 4 
		merged_seq, merged_qual = ok_overlap(loffset, overl, fseq, fqual, rseq, rqual)
	elif len(fseq) - int(far_report[QEND]) > max_align_offset or int(far_report[SSTART]) > max_align_offset:
		"""(randomly occured alignment in center of sequences, need to blast)
		"""
		# faref means forward [read] against reference [sequence]
		# raref means reverse [read] against reference [sequence]
		gap = False
		for faref_report, raref_report in align_reads_ag_some_ref(rseq):
			if (faref_report, raref_report) == ("NOTHING", "MORE"):
				break
			forw_start = int(faref_report[SSTART]) - int(faref_report[QSTART])
			forw_end = forw_start + len(fseq)
			rev_start = int(raref_report[SSTART]) - int(raref_report[QSTART])
			rev_end = rev_start + len(rseq)

			if forw_end < rev_start:
				gap = True
				break

		if not gap:
			print("\nPUTATIVE CHIMERA!!!")
			# WRITE FASTQ RECORD
		else:
			print("\nGAP")
			# fill the gap with N
			gap_len = rev_start - forw_end
			merged_seq = fseq + 'N' * (gap_len - 1) + rseq
			merged_qual = fqual + chr(33) * (gap_len - 1) + rqual
			print(merged_seq)
	else:
		print("\nWeird reads!")
		"""(putative artifact)
		  FFFFFFF
		RRRRRRR
		"""
		# WRITE FASTQ RECORD



	# merged_file.write(fastq_recs["mR1"]["seq_id"] + '\n')
	# merged_file.write(merged_seq + '\n')
	# merged_file.write(fastq_recs["mR1"]["optional_id"] + '\n')
	# merged_file.write(merged_qual + '\n')


for key in filt_read_files.keys():
	filt_read_files[key].close()