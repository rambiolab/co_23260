import glob

SAMPLES = ["ERR4398741", "ERR4398768", "ERR4398769"]


rule all:
	input:
		expand("results/fastqc/{sample}_1_fastqc.html", sample=SAMPLES),
		expand("results/fastqc/{sample}_2_fastqc.html", sample=SAMPLES),
		expand("results/bbduk/{sample}_1.bbduk.trim.fastq.gz", sample=SAMPLES),
		expand("results/bbduk/{sample}_2.bbduk.trim.fastq.gz", sample=SAMPLES),
		expand("results/bbduk/{sample}_s.bbduk.trim.fastq.gz", sample=SAMPLES),
		expand("results/fastp/{sample}_1.fastp.trim.fastq.gz", sample=SAMPLES),
		expand("results/fastp/{sample}_2.fastp.trim.fastq.gz", sample=SAMPLES),
		expand("results/fastp/{sample}_s.fastp.trim.fastq.gz", sample=SAMPLES),
		expand("results/fastp/{sample}.html", sample=SAMPLES),
		expand("results/kraken/{sample}_kraken.out", sample=SAMPLES),
		expand("results/kraken2/{sample}_kraken2.out", sample=SAMPLES),
		expand("results/kraken2/{sample}_kraken2.report", sample=SAMPLES),
		expand("results/bracken/{sample}_bracken.S.out", sample=SAMPLES),
		expand("results/bracken/{sample}_bracken.S.report", sample=SAMPLES),
		expand("results/bracken/{sample}_bracken.S.summary", sample=SAMPLES),
		expand("results/kma_resfinder/{sample}/{sample}.mapstat", sample=SAMPLES),
		expand("results/kma_silva/{sample}/{sample}.mapstat", sample=SAMPLES),
		expand("results/assembly/{sample}/scaffolds.fasta", sample=SAMPLES),
		expand("results/assembly/{sample}/{sample}_scaffolds_filtered.fasta", sample=SAMPLES),
		expand("results/mapped/{sample}/mapped.sam",sample=SAMPLES),
		expand("results/mapped/{sample}/mapped.sort.bam",sample=SAMPLES),
		expand("results/mapped/{sample}/coverage.txt", sample=SAMPLES),
		"results/checkm/quality_report.tsv",
		"results/checkm/all_bin_list.txt",
		"results/checkm/edited_quality_report.txt",
		"results/drep/drep.done", 
		"results/drep/drep2.done",
		"results/gtdbtk/classify.done"	

rule bbduk:
	envmodules:
		"tools",
		"jdk/19",
		"bbmap/38.90"
	input:
		infile1="raw_reads/{sample}_1.fastq.gz",
		infile2="raw_reads/{sample}_2.fastq.gz"
	output:
		outfile1="results/bbduk/{sample}_1.bbduk.trim.fastq.gz",
		outfile2="results/bbduk/{sample}_2.bbduk.trim.fastq.gz",
		singletons="results/bbduk/{sample}_s.bbduk.trim.fastq.gz"
	params:
		minlength=50,
		time="results/bbduk/time.{sample}.log"
	threads:
		10
	resources:
		mem_gb=40
	shell:
		"""
		/usr/bin/time -v -o {params.time} bbduk.sh in={input.infile1} in2={input.infile2} out={output.outfile1} out2={output.outfile2} outs={output.singletons} threads={threads} minlength={params.minlength} k=19 kmin=11 ktrim=r qtrim=r trimq=20 minlength=50 overwrite=t tbo ziplevel=6 ref=/home/projects/co_23260/data/adaptors/adaptors.fa
		"""

rule fastp:
	envmodules:
		"tools",
		"fastp/0.23.2"
	input:
		infile1="raw_reads/{sample}_1.fastq.gz",
		infile2="raw_reads/{sample}_2.fastq.gz"
	output:
		out="results/fastp/{sample}_1.fastp.trim.fastq.gz",
		out2="results/fastp/{sample}_2.fastp.trim.fastq.gz",
		outs="results/fastp/{sample}_s.fastp.trim.fastq.gz",
		html="results/fastp/{sample}.html",
		json="results/fastp/{sample}.json"
	params:
		minlength=50,
		outm="results/fastp/{sample}_m.fastp.trim.fastq.gz",
		time="results/fastp/time.{sample}.log"
	threads:
		10
	resources:
		mem_gb=20
	shell:
		"""
		/usr/bin/time -v -o {params.time} fastp -i {input.infile1} -I {input.infile2} -o {output.out} -O {output.out2} --merge --merged_out {params.outm} --unpaired1 {output.outs} --unpaired2 {output.outs}  --thread {threads} --average_qual 20 --length_required {params.minlength} --cut_tail -h {output.html} -j {output.json} --overlap_diff_limit 1 --detect_adapter_for_pe
		cat {params.outm} >> {output.outs}
		rm {params.outm}
		"""

rule fastqc:
	envmodules:
		"tools",
		"perl/5.30.2",
		"jdk/19",
		"fastqc/0.11.9"
	input:
		"raw_reads/{sample}_1.fastq.gz",
		"raw_reads/{sample}_2.fastq.gz",
		"results/bbduk/{sample}_1.bbduk.trim.fastq.gz",
		"results/bbduk/{sample}_2.bbduk.trim.fastq.gz",
		"results/bbduk/{sample}_s.bbduk.trim.fastq.gz",
		"results/fastp/{sample}_1.fastp.trim.fastq.gz",
		"results/fastp/{sample}_2.fastp.trim.fastq.gz",
		"results/fastp/{sample}_s.fastp.trim.fastq.gz"
	output:
		"results/fastqc/{sample}_1_fastqc.html",
		"results/fastqc/{sample}_2_fastqc.html"
	params:
		outdir="results/fastqc",
		time="results/fastqc/time.{sample}.log"
	threads:
		10
	resources:
		mem_gb=10
	shell:
		"""
		/usr/bin/time -v -o {params.time} fastqc --threads {threads} -o {params.outdir} {input}
		"""

rule kraken:
	envmodules:
		"tools",
		"kraken/1.1"
	input:
		in1="results/bbduk/{sample}_1.bbduk.trim.fastq.gz",
		in2="results/bbduk/{sample}_2.bbduk.trim.fastq.gz"
	output:
		out="results/kraken/{sample}_kraken.out",
		report="results/kraken/{sample}_kraken.report"
	params:
		time="results/kraken/time.{sample}.log",
		db="/home/projects/co_23260/data/databases/minikraken_20171013_4GB"
	threads:
		8
	resources:
		mem_gb=30
	shell:
		"""
		/usr/bin/time -v -o {params.time} kraken --gzip-compressed --fastq-input --threads {threads} --db {params.db} --output {output.out} --paired {input.in1} {input.in2}
		kraken-report -db {params.db} {output.out} > {output.report}
		"""

rule kraken2:
	envmodules:
		"tools",
		"kraken2/20231104"
	input:
		in1="results/bbduk/{sample}_1.bbduk.trim.fastq.gz",
		in2="results/bbduk/{sample}_2.bbduk.trim.fastq.gz"
	output:
		out="results/kraken2/{sample}_kraken2.out",
		report="results/kraken2/{sample}_kraken2.report"
	params:
		time="results/kraken2/time.{sample}.log",
		db="/home/projects/co_23260/data/databases/minikraken2_v2_8GB_201904_UPDATE"
	threads:
		8
	resources:
		mem_gb=62
	shell:
		"""
		/usr/bin/time -v -o {params.time} kraken2 --report {output.report} --gzip-compressed --threads {threads} --db {params.db} --output {output.out} --paired {input.in1} {input.in2}
		"""

rule bracken:
	envmodules:
		"tools",
		"bracken/3.1"
	input:
		inp="results/kraken2/{sample}_kraken2.report"
	output:
		out="results/bracken/{sample}_bracken.S.out",
		report="results/bracken/{sample}_bracken.S.report",
		summary="results/bracken/{sample}_bracken.S.summary"
	params:
		time="results/bracken/time.{sample}.log",
		db="/home/projects/co_23260/data/databases/minikraken2_v2_8GB_201904_UPDATE"
	shell:
		"""
		/usr/bin/time -v -o {params.time} bracken -d {params.db} -i {input.inp} -o {output.out} -w {output.report} -r 150 -l S -t 10 > {output.summary}
		"""

rule kma_resfinder:
	input:
		R1="results/bbduk/{sample}_1.bbduk.trim.fastq.gz",
		R2="results/bbduk/{sample}_2.bbduk.trim.fastq.gz",
		S="results/bbduk/{sample}_s.bbduk.trim.fastq.gz"
	output:
		"results/kma_resfinder/{sample}/{sample}.mapstat"
	params:
		outdir="results/kma_resfinder/{sample}/{sample}",
		time="results/kma_resfinder/time.{sample}.log",
		db="/home/projects/co_23260/data/groups/tnpe/resfinder_db/ResFinder",
		threads=1,
		kma="/home/projects/co_23260/data/groups/tnpe/kma/kma"
	resources:
		mem_gb=10
	shell:
		"""
		/usr/bin/time -v -o {params.time} {params.kma} -ipe {input.R1} {input.R2} -i {input.S} -t_db {params.db} -o {params.outdir}  -ef -reassign  -ill -t {params.threads}
		"""

rule kma_silva:
	input:
		R1="results/bbduk/{sample}_1.bbduk.trim.fastq.gz",
		R2="results/bbduk/{sample}_2.bbduk.trim.fastq.gz",
		S="results/bbduk/{sample}_s.bbduk.trim.fastq.gz"
	output:
		"results/kma_silva/{sample}/{sample}.mapstat"
	params:
		outdir="results/kma_silva/{sample}/{sample}",
		time="results/kma_silva/time.{sample}.log",
		db="/home/projects/co_23260/data/databases/kma_silva_v138/Silva_20200116",
		threads=12,
		kma="/home/projects/co_23260/data/groups/tnpe/kma/kma"
	resources:
		mem_gb=60
	shell:
		"""
		/usr/bin/time -v -o {params.time} {params.kma} -ipe {input.R1} {input.R2} -i {input.S} -t_db {params.db} -o {params.outdir}  -ef -reassign  -ill -t {params.threads} -mem_mode
		"""

rule assembly:
	envmodules:
		"tools",
		"anaconda3/2025.06-1",
		"spades/4.0.0"
	input:
		in1="results/bbduk/{sample}_1.bbduk.trim.fastq.gz",
		in2="results/bbduk/{sample}_2.bbduk.trim.fastq.gz",
		singletons="results/bbduk/{sample}_s.bbduk.trim.fastq.gz"
	output:
		"results/assembly/{sample}/scaffolds.fasta"
	params:
		outdir="results/assembly/{sample}",
		time="results/assembly/{sample}/time.log"
	threads:
		20
	resources:
		mem_gb=180
	shell:
		"""
		/usr/bin/time -v -o {params.time} metaspades.py -t {threads} -1 {input.in1} -2 {input.in2} -s {input.singletons}  -k 27,47,67,87,107,127 --memory {resources.mem_gb} -o {params.outdir}
		"""

rule filter:
	envmodules:
		"tools",
		"jdk/19",
		"bbmap/38.90"
	input:
		infile="results/assembly/{sample}/scaffolds.fasta"
	output:
		outfile="results/assembly/{sample}/{sample}_scaffolds_filtered.fasta"
	params:
		pfx="{sample}"
	shell:
		"""
		rename.sh ow=t fastawrap=60 minscaf=1000 prefix={params.pfx} addprefix=t in={input.infile} out={output.outfile}
		"""

rule bbmap:
	envmodules:
		"tools",
		"jdk/19",
		"bbmap/38.90"
	input:
		infile="results/bbduk/{sample}_1.bbduk.trim.fastq.gz",
		infile2="results/bbduk/{sample}_2.bbduk.trim.fastq.gz", 
		ref="results/assembly/{sample}/{sample}_scaffolds_filtered.fasta"
	output:
		sam="results/mapped/{sample}/mapped.sam"
	threads:
		40
	resources:
		mem_gb=180
	shell:
		"""
		bbmap.sh in={input.infile} in2={input.infile2} minid=0.97 threads={threads} ref={input.ref} outm={output.sam} overwrite=t nodisk=t

		"""


rule samtools:
	envmodules:
		"tools",
		"samtools/1.20"
	input:
		infile = "results/mapped/{sample}/mapped.sam"
	output:
		outfile ="results/mapped/{sample}/mapped.sort.bam"
	shell:
		"""
		samtools view -bSh1 {input.infile} | samtools sort -m 20G -@ 3 > {output.outfile}
               
		"""


rule jgi:
	envmodules:
		"ngs",
		"tools",
		"perl/5.24.0",
		"metabat/2.17",
		"samtools/1.20",
		"java/17-openjdk",
		"bbmap/38.90"
	input:
		infile = "results/mapped/{sample}/mapped.sort.bam"
	output:
		outfile = "results/mapped/{sample}/coverage.txt"
	shell:
		"""
		jgi_summarize_bam_contig_depths {input.infile} --outputDepth {output.outfile}

		"""

rule metabat2:
	envmodules:
		"tools",
		"perl/5.24.0",
		"metabat/2.12.1"		
	input:
		infile1="results/assembly/{sample}/{sample}_scaffolds_filtered.fasta",
		infile2="results/mapped/{sample}/coverage.txt"
	output:
		outdir = directory("results/metabat2/{sample}")
	threads:
		20
	shell:
		"""
		mkdir -p {output.outdir}
		metabat2 -i {input.infile1} -a {input.infile2} -o {output.outdir}/{wildcards.sample}.bin -t {threads}
		"""

rule checkm2:
	envmodules: 
		"tools",
		"anaconda3/2023.09-0",
		"checkm2/1.0.2"
	input:
		indir=expand("results/metabat2/{sample}", sample=SAMPLES)
	output:
		outfile="results/checkm/quality_report.tsv"
	params:
		outdir="results/checkm",
		indir="results/all_bins",
		database="/home/projects/co_23260/data/databases/CheckM2_database/uniref100.KO.1.dmnd"
	threads: 
		40
	shell:
		"""
		mkdir -p {params.indir}
                ln -sf $(pwd)/results/metabat2/*/*.fa {params.indir}
		checkm2 predict  -x fa -i {params.indir} -o {params.outdir} -t {threads} --database_path {params.database}
		"""


rule make_drep_fofn:
	input:
		flag="results/checkm/quality_report.tsv"
	output:
		outfile="results/checkm/all_bin_list.txt"
	shell:
		"""
		mkdir -p results/drep
		ls -1 $(pwd)/results/all_bins/*.fa > {output.outfile}
		"""

rule drep_quality_table: 
	envmodules: 
		"tools",
		"anaconda3/2025.06-1",
		"mash/2.3",
		"fastani/1.34",
		"prodigal/2.6.3"
	input:
		infile = "results/checkm/quality_report.tsv"
	output:
		outfile="results/checkm/edited_quality_report.txt"
	params:
		script="/home/projects/co_23260/apps/edit_checkm2_results.py"
	shell:
		"""
		python3 {params.script} {input.infile}
		mv edited_quality_report.txt {output.outfile}
		"""


rule drep: 
	envmodules:
		"tools",
		"anaconda3/2025.06-1",
		"mash/2.3",
		"fastani/1.34",
		"prodigal/2.6.3" 
	input: 
		infile1 = "results/checkm/all_bin_list.txt",
		infile2 = "results/checkm/edited_quality_report.txt"
	output: 
		outdir = directory("results/drep"),
		flag = "results/drep/drep.done"
	threads: 
		40	
	shell: 
		"""	
		mkdir -p {output.outdir}
		dRep dereplicate {output.outdir} -g {input.infile1} -p {threads} -pa 0.9 -sa 0.95 --genomeInfo {input.infile2}
		touch {output.flag}
		"""





rule dRep_another_way:
	envmodules:
		"ngs tools",
		"anaconda3/2025.06-1",
		"mash/2.3",
		"prodigal/2.6.3",
		"centrifuge/1.0.4-beta",
		"hmmer/3.1b2",
		"pplacer/1.1.alpha17",
		"checkm2/20220719",
		"mummer/4.0.1",
		"fastani/1.34"
	input: 
		infile1 = "results/checkm/all_bin_list.txt",
		infile2 = "results/checkm/edited_quality_report.txt"
	output:
		outdir = directory("results/drep2"),
		flag = "results/drep/drep2.done"
	threads: 
		40
	shell:
		"""
		mkdir -p {output.outdir}
		dRep dereplicate {output.outdir} -g {input.infile1} -p {threads} -pa 0.9 -sa 0.95 --genomeInfo {input.infile2}
		touch {output.flag}
		"""

rule gtdb_tk:
	envmodules: 
		"tools",
		"gtdbtk/2.1.1"
	input:
		infile="results/drep/drep2.done"
	output: 
		flag= "results/gtdbtk/classify.done",
	params:
		outdir= "results/gtdbtk", 
		indir="results/drep2/dereplicated_genomes/"
	threads: 
		40
	shell:		
		"""
		gtdbtk classify_wf --genome_dir {params.indir} --out_dir {params.outdir} --extension fa --cpus {threads}
		touch {output.flag}
		"""
	
