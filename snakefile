## =============================================================================
## WORKFLOW PROJECT: RRBS
## INIT DATE: 2023
import pandas as pd
import os, os.path, sys
from os.path import join, abspath
from snakemake.utils import validate#, min_version

## If SettingwithCopyWarning is ever generated, raise error instead of warning; makes code more strict
pd.set_option("mode.chained_assignment", "raise")

## DEFAULT PATHS
DIR_ENVS    = abspath("envs")
DIR_SCHEMAS = abspath("schemas")

## LOAD VARIABLES FROM CONFIGFILE
## submit on command-line via --configfile
if config=={}:
    sys.stderr.write("Please submit config-file with --configfile <file>. Exit.\n")
    sys.exit(1)

bs_sum_file = config["bs_summary"]["bs_summary_file"]

trim_params = None
if config["trimming"]["perform"]:
    trim_params = config["trimming"]["extra"]

phi_index = None
if config["phi_rem"]["perform"]:
    phi_index = config["phi_rem"]["index"]

#print(trim_params)
#print(type(trim_params))
#print(phi_index)

## Setup result dirs
DIR_BASE       = config["resultdir"]
DIR_LOGS       = join(DIR_BASE, "logs")
DIR_BENCHMARKS = join(DIR_BASE, "benchmarks")
DIR_RES        = join(DIR_BASE, "results")

SAMPLES = pd.read_csv(abspath(config["samples"]), sep="\t", keep_default_na=False, na_values="").set_index("sample", drop=False)
validate(SAMPLES, schema=join(DIR_SCHEMAS, "samples.schema.yaml"))

## reading samplename from samplesheet
sys.stderr.write("Reading samples from samplesheet: '{}' ...\n".format(config["samples"]))

## test if sample configuration is correct (expected paired-end)
for fname in SAMPLES["sample"]:

    fname1 = pd.notnull(SAMPLES.loc[fname, "read1"])  # True or False
    fname2 = pd.notnull(SAMPLES.loc[fname, "read2"]) # True or False

    if fname1:
        if not os.path.isfile(SAMPLES.loc[fname, "read1"]):
            sys.stderr.write("File '{}: read1' from samplesheet can not be found. Make sure the file exists. Exit\n".format(fname))
            sys.exit(1)
        if not fname2:
            sys.stderr.write("For {} PE, please specify path to forward and reverse reads in columns 1 and 2, respectively\n".format(fname))
            sys.exit(1)
        else:
            if not os.path.isfile(SAMPLES.loc[fname, "read2"]):
                sys.stderr.write("File '{}: read2' from samplesheet can not be found. Make sure the file exists. Exit\n".format(fname))
                sys.exit(1)
    else:
        sys.stderr.write("For {} PE, please specify path to forward and reverse reads in columns 1 and 2, respectively\n".format(fname))
        sys.exit(1)

#print(os.path.dirname("./"))
## Combinations of Input Read Files Passed Checks     
sys.stderr.write("Combinations of Input Read Files Passed Checks\n")

NUM_SAMPLES = len(SAMPLES["sample"])
sys.stderr.write("{} samples to process\n".format(NUM_SAMPLES))
sys.stderr.flush()

## =============================================================================
## FUNCTIONS FOR RULE INPUTS
## =============================================================================
def get_fq_pe(wildcards):
    return [SAMPLES.loc[wildcards.sample, "read1"], 
            SAMPLES.loc[wildcards.sample, "read2"]]

def get_data_for_phiX(wildcards):
    if config["trimming"]["perform"]:
        return [join(DIR_RES, "trimmed/{0}/{0}_val_1.fq.gz".format(wildcards.sample)),
                join(DIR_RES, "trimmed/{0}/{0}_val_2.fq.gz".format(wildcards.sample))]
    else:
        return get_fq_pe(wildcards)


def get_data_for_bismark(wildcards):
    if config["phi_rem"]["perform"]:
        # phiX removed
        return [join(DIR_RES, "phix_rem/{0}/{0}_phixrem.1.fq.gz".format(wildcards.sample)),
                join(DIR_RES, "phix_rem/{0}/{0}_phixrem.2.fq.gz".format(wildcards.sample))]
    elif config["trimming"]["perform"]:
#        # trimmed
        return [join(DIR_RES, "trimmed/{0}/{0}_val_1.fq.gz".format(wildcards.sample)),
                join(DIR_RES, "trimmed/{0}/{0}_val_2.fq.gz".format(wildcards.sample))]
    else:
         # raw 
        return get_fq_pe(wildcards)

bisfiles = expand(join(DIR_RES, "bismark/{sample}/{sample}_bs_pe.bam"), sample=SAMPLES["sample"])
## =============================================================================
## RULES
## =============================================================================

#TARGETS = expand([join(DIR_RES, "{sample}.txt"), join(DIR_RES, "{sample}.done")], sample=SAMPLES["sample"])

rule all:
    input: 
        expand(join(DIR_RES, "bismark/{sample}/meth_extract/{sample}_bs_pe.bismark.cov.gz"),  sample=SAMPLES["sample"]), 
        join(join(DIR_RES, "bismark"), "{}.txt".format(bs_sum_file)), join(join(DIR_RES, "bismark"), "{}.html".format(bs_sum_file))
    

## 1. trim with trim galore which has rrbs capabilities. note that --gzip is added to make sure that no matter the input file, output is gzipped for pipeline compatibility
## note for output: in a forum, they say it is okay even if the output is not explicitly in rule's command; snakemake will look for input-output chaining
rule rrbs_trim:
    input:
        fq = get_fq_pe
    output:
        trim_dir = directory(join(DIR_RES, "trimmed/{sample}")),
        trim1 = temp(join(DIR_RES, "trimmed/{sample}/{sample}_val_1.fq.gz")),
        trim2 = temp(join(DIR_RES, "trimmed/{sample}/{sample}_val_2.fq.gz"))
    priority: 10
    log:
        join(DIR_LOGS, "trimmed/{sample}_trim_pe.log")
    benchmark:
        join(DIR_BENCHMARKS, "trimmed/{sample}_trim_pe.txt")
    threads: 4
    params:
        trim_params = trim_params,
        bn = "{sample}"
    conda:
        join(DIR_ENVS, "rrbs.yaml")
    shell:
        "nice trim_galore --cores {threads} {params.trim_params} --paired --fastqc --gzip --basename {params.bn} {input.fq[0]} {input.fq[1]} -o {output.trim_dir} > {log} 2>&1"

## 2. remove phix sequences

rule phix_rem:
    input:
        btin = get_data_for_phiX
    output:
        btsam = temp(join(DIR_RES, "phix_rem/{sample}/{sample}.sam")),
        btun1 = join(DIR_RES, "phix_rem/{sample}/{sample}_phixrem.1"),
	    btun2 = join(DIR_RES, "phix_rem/{sample}/{sample}_phixrem.2")
    priority: 20
    log:
        join(DIR_LOGS, "phix_rem/{sample}_phixrem.log")
    benchmark:
        join(DIR_BENCHMARKS, "phix_rem/{sample}_phixrem.txt")
    threads: 4
    params:
        phi_index=phi_index,
	btun = join(DIR_RES, "phix_rem/{sample}/{sample}_phixrem")
    conda:
        join(DIR_ENVS, "rrbs.yaml")
    shell:
        "nice bowtie2 -p {threads} -x {params.phi_index} -1 {input.btin[0]} -2 {input.btin[1]} --un-conc-gz {params.btun} -S {output.btsam} > {log} 2>&1"

## 2.1 ##phix unmapped outputs doesn't come in gz, so just rename to have .gz
rule rename_phi_un:
    input: 
        rules.phix_rem.output.btun1,
	rules.phix_rem.output.btun2
    priority: 20
    output:
        un1 = temp(join(DIR_RES, "phix_rem/{sample}/{sample}_phixrem.1.fq.gz")),
        un2 = temp(join(DIR_RES, "phix_rem/{sample}/{sample}_phixrem.2.fq.gz"))
    shell:
        """
        mv {input[0]}  {output.un1}
        mv {input[1]}  {output.un2}
        """

## 3. bismark run ; constrain resources because they said this could blow up cores and I don't want to blow servers up when they are used simultaneously

rule bismark_main:
    input: 
        bisin = get_data_for_bismark
    output: 
        odir = directory(join(DIR_RES, "bismark/{sample}")),
        obam = temp(join(DIR_RES, "bismark/{sample}/{sample}_bs_pe.bam"))
    priority: 30
    log:
        join(DIR_LOGS, "bismark/{sample}_bismark.log")
    benchmark:
        join(DIR_BENCHMARKS, "bismark/{sample}_bismark.txt")
    params:
        bs_index = config["bs_main"]["bs_index"],
        bs_extra = config["bs_main"]["bs_extra"],
        bs_basename = "{sample}"
    conda:
        join(DIR_ENVS, "rrbs.yaml")
    resources:
        const=1
    shell: 
        "nice bismark --genome {params.bs_index} -1 {input.bisin[0]} -2 {input.bisin[1]} {params.bs_extra} -B {params.bs_basename}_bs --temp_dir {output.odir} -o {output.odir} > {log} 2>&1"

# 4. bismark methylation extraction; constrain resources because they said this could blow up cores and I don't want to blow servers up when they are used simultaneously


rule bismark_meth:
    input:
        rules.bismark_main.output.obam
    output: 
        bsmethdir = directory(join(DIR_RES, "bismark/{sample}/meth_extract")),
        bdcovfiles = join(DIR_RES, "bismark/{sample}/meth_extract/{sample}_bs_pe.bismark.cov.gz")
    #priority: 30
    log: 
        join(DIR_LOGS, "bismark/{sample}_bismark_meth.log")
    benchmark:
        join(DIR_BENCHMARKS, "bismark/{sample}_bismark_meth.txt")
    params: 
        bs_meth_extra = config["bs_meth"]["bs_meth_extra"]
    conda:
        join(DIR_ENVS, "rrbs.yaml")
    resources:
        const=1
    shell: 
        "nice bismark_methylation_extractor {params.bs_meth_extra} --bedGraph -p -o {output.bsmethdir} {input} > {log} 2>&1"

rule bismark_summary:
    input: 
        bisfiles
    output:
        join(join(DIR_RES, "bismark"), "{}.txt".format(bs_sum_file)),
        join(join(DIR_RES, "bismark"), "{}.html".format(bs_sum_file))
    log:
        join(DIR_LOGS, "bismark/summaryfile.log")
    benchmark:
        join(DIR_BENCHMARKS, "bismark/summaryfile.txt")
    params:
        join(join(DIR_RES, "bismark"), bs_sum_file)
    conda:
        join(DIR_ENVS, "rrbs.yaml")
    shell:
        "nice bismark2summary {input} -o {params} > {log} 2>&1"