


# RNA-Seq Reads to Counts Translation

## Introduction

This tutorial demonstrates translation of a Galaxy workflow to Nextflow using `janis translate`. 

<br>


**Source Workflow**

The workflow used in this tutorial is taken from the [McDonnell Genome Institute](https://www.genome.wustl.edu/) (MGI) [analysis-workflows](https://github.com/genome/analysis-workflows) repository. <br>
This resource stores publically available analysis pipelines for genomics data. <br>
It is a fantastic piece of research software, and the authors thank MGI for their contribution to open-source research software. 

The workflow using in this tutorial - [align_sort_markdup](https://github.com/genome/analysis-workflows/blob/master/definitions/subworkflows/align_sort_markdup.cwl) - accepts multiple unaligned readsets as input and produces a single polished alignment bam file. 

*Main Inputs* 
- Unaligned reads stored across multiple BAM files
- Reference genome index

*Main Outputs*
- Single BAM file storing alignment for all readsets

*Steps*
- Read alignment (run in parallel across all readsets) - `bwa mem`
- Merging alignment BAM files to single file - `samtools merge`
- Sorting merged BAM by coordinate - `sambamba sort`
- Tagging duplicate reads in alignment - `picard MarkDuplicates` 
- Indexing final BAM - `samtools index`

<br>

**Tutorial Outcomes**

In this tutorial we will:
- Install the required software
- Translate the CWL using `janis translate`
- Make manual adjustments to the translation if necessary
- Run the nextflow using sample input data to validate our nextflow code

After completing this short tutorial, you will be familiar with using `janis translate` to migrate workflow tools in CWL to Nextflow.

Other tutorials exist to demonstrate migration from WDL / CWL / Galaxy -> Nextflow in this repository, including full workflow migrations with multiple tasks. 

<br>

**Installation**

To begin, make sure you have [nextflow](https://nf-co.re/usage/installation), [docker](https://docs.docker.com/engine/install/), and [janis translate](https://janis.readthedocs.io/en/latest/index.html) installed. <br>
The links above contain installation instructions. 

<br>

## Janis Translate

To translate a workflow,  we use `janis translate`.

```
janis translate --from <src> --to <dest> <filepath>
```

The `--from` specifies the workflow language of the source file(s), and `--to` specifies the destination we want to translate to. 

<br>

**Run Janis Translate**

In our case, we want to translate CWL -> Nextflow, and our source CWL file is located at `source/subworkflows/align_sort_markdup.cwl` relative to this document.

*using pip*

To translate `align_sort_markdup.cwl` to nextflow, we can write the following in a shell:
```
janis translate --from cwl --to nextflow ./source/subworkflows/align_sort_markdup.cwl
```

*using docker (linux bash)*

If the janis translate docker container is being used, we can write the following:
```
docker run -v $(pwd):/home janis translate --from cwl --to nextflow ./source/subworkflows/align_sort_markdup.cwl
```

<br>

**Translation Output**

The output translation will contain multiple files and directories.<br>
You will see a folder called `translated` appear - inside this folder, we should see the following structure:

```
├── main.nf                             # main workflow (align_sort_markdup)
├── modules                             # folder containing nextflow processes
│   ├── align_and_tag.nf
│   ├── index_bam.nf
│   ├── mark_duplicates_and_sort.nf
│   ├── merge_bams_samtools.nf
│   └── name_sort.nf
├── nextflow.config                     # config file to supply input information
├── subworkflows                        # folder containing nextflow subworkflows
│   └── align.nf    
└── templates                           # folder containing any scripts used by processes
    └── markduplicates_helper.sh
```

Now we have performed translation using `janis translate`, we need to check the translated workflow for correctness.  

From here, we will do a test-run of the workflow using sample data, and make manual adjustments to the translated workflow where needed. 

<br>

## Running the Translated Workflow


**Inspect main.nf**

The main workflow translation appears as `main.nf` in the `translated/` folder. <br>

This filename is just a convention, and we use it to provide clarity about the main entry point of the workflow. <br>
In our case `main.nf` holds the nextflow definition for the  `align_sort_markdup.cwl` workflow. 

> NOTE: <br>
> Before continuing, feel free to have a look at the other nextflow files which have been generated during translation:<br>
> Each CWL subworkflow appears as a nextflow `workflow` in the `subworkflows/` directory.<br>
> Each CWL tool appears as a nextflow `process` in the `modules/` directory. 

In `main.nf` we see the nextflow workflows / processes called by the main workflow:

```
include { ALIGN } from './subworkflows/align'
include { INDEX_BAM } from './modules/index_bam'
include { MARK_DUPLICATES_AND_SORT } from './modules/mark_duplicates_and_sort'
include { MERGE_BAMS_SAMTOOLS as MERGE } from './modules/merge_bams_samtools'
include { NAME_SORT } from './modules/name_sort'
```

We also see that some nextflow `Channels` have been set up. <br>
These are used to supply data according to nextflow's adoption of the *dataflow* programming model.
```
// data which will be passed as channels
ch_bams        = Channel.fromPath( params.bams ).toList()
ch_reference   = Channel.fromPath( params.reference ).toList()
ch_readgroups  = Channel.of( params.readgroups ).toList()
```

Focusing on the channel declarations, we want to note a few things:

- `ch_bams` is analygous to the *'bams'* input in `align_sort_markdup.cwl`. <br>
It declares a queue channel which expects the data supplied via `params.bams` are `path` types. <br>
It then groups the bams together as a sole emission. <br>
We will need to set up `params.bams` to supply this data.

- `ch_reference` is analygous to the *'reference'* input in `align_sort_markdup.cwl`. <br>
This channel collects the *'reference'* primary & secondary files as a sole emission. <br> 
We will need to set up `params.reference` to supply this data. 

- `ch_readgroups` is analygous to the *'readgroups'* input in `align_sort_markdup.cwl`. <br>
It is the same as `ch_bams`, except it requires `val` types rather than `path` types. <br>
We will need to set up `params.readgroups` to supply this data.


> Note: `.toList()` <br><br>
> Nextflow queue channels work differently to lists. <br>
> Instead of supplying all items together, queue channels emit each item separately. <br> 
> This results in a separate task being spawned for each item in the queue when the channel is used. <br>
> As the CWL workflow input specifies that `bams` is a list, we use `.toList()` to group all items as a sole emission. <br>
> This mimics a CWL array which is the datatype of the `bams` inputs. <br><br>
> As it turns out, the CWL workflow ends up running the `align` step in parallel across the `bams` & `readgroups` inputs. <br><br>
> Parallelisation in nextflow happens by default. <br>
> To facilitate this, the `.flatten()` method is called on  `ch_bams` and `ch_readgroups` when used in the `ALIGN` task. <br>
> This emits items in `ch_bams` and `ch_readgroups` individually, spawning a new `ALIGN` task for each pair. <br><br>
> We're kinda doing redundant work by calling `.toList()`, then `.flatten()` when `ch_bams` and `ch_readgroups` are used.<br>
> `janis translate` isn't smart enough yet to detect this yet, but may do so in future. 

<br>

The main `workflow {}` has 5 tasks. <br>
Each task has been supplied values according to the source workflow. <br>
Comments display the name of the process/workflow input which is being fed a particular value. 

```
workflow {

    ALIGN(
        ch_bams.flatten(),       // bam
        ch_reference,            // reference
        ch_readgroups.flatten()  // readgroup
    )

    INDEX_BAM(
        MARK_DUPLICATES_AND_SORT.out.sorted_bam.map{ tuple -> tuple[0] }  // bam
    )

    MARK_DUPLICATES_AND_SORT(
        params.mark_duplicates_and_sort.script,  // script
        NAME_SORT.out.name_sorted_bam,           // bam
        params.NULL_VALUE,                       // input_sort_order
        params.final_name                        // output_name
    )

    MERGE(
        ALIGN.out.tagged_bam.toList(),  // bams
        params.final_name               // name
    )

    NAME_SORT(
        MERGE.out.merged_bam  // bam
    )

}
```

<br>

Before `main.nf` can be run, we will need to supply values for the `params` variables. 
This is done in `nextflow.config`. 

<br>

**Inspect nextflow.config**

To test the translated workflow, we will first set up workflow inputs in `nextflow.config`. 

Before running a workflow, nextflow will attempt to open `nextflow.config` and read in config information and global *param* variables from this file. 
We use this file to tell nextflow how to run and to supply workflow inputs.

Inside the `translated/` folder you will see that `nextflow.config` is already provided. 

Janis translate creates this file to provide clarity about the necessary workflow inputs, and to set some other config variables. 

Open `nextflow.config` and have a look at the contents. It should look similar to the following: 

```
nextflow.enable.dsl=2
docker.enabled = true

params {
    
    // Placeholder for null values.
    // Do not alter unless you know what you are doing.
    NULL_VALUE = 'NULL'

    // WORKFLOW OUTPUT DIRECTORY
    outdir  = './outputs'

    // INPUTS (MANDATORY)
    bams        = []          // (MANDATORY array)             eg. [file1, ...]
    reference   = []          // (MANDATORY fastawithindexes)  eg. [fasta, amb, ann, bwt, dict, fai, pac, sa]
    readgroups  = NULL_VALUE  // (MANDATORY array)             eg. [string1, ...]

    // INPUTS (OPTIONAL)
    final_name  = "final.bam" 

    // PROCESS: ALIGN_AND_TAG
    align_and_tag.cpus    = 8     
    align_and_tag.memory  = 20000 

    // PROCESS: INDEX_BAM
    index_bam.memory  = 4000 

    // PROCESS: MARK_DUPLICATES_AND_SORT
    mark_duplicates_and_sort.script  = "(local_dir)/templates/markduplicates_helper.sh"
    mark_duplicates_and_sort.cpus    = 8
    mark_duplicates_and_sort.memory  = 40000

    // PROCESS: MERGE_BAMS_SAMTOOLS
    merge_bams_samtools.cpus    = 4    
    merge_bams_samtools.memory  = 8000 

    // PROCESS: NAME_SORT
    name_sort.cpus    = 8     
    name_sort.memory  = 26000 

}
```

> NOTE: `NULL_VALUE = 'NULL'`<br><br>
> Nextflow doesn't like `null` values to be passed to process inputs. <br>
> This is a challenge for translation as other languages allow `optional` inputs. <br>
> To get around this, Janis Translate sets the `params.NULL_VALUE` variable as a `null` placeholder for `val` type inputs. <br>
> You will see this being used in nextflow processes to do optionality checking.

<br>

The auto-generated `nextflow.config` splits up workflow inputs using some headings. 

```
// INPUTS (MANDATORY)
```
Workflow inputs which are required to run the workflow. We must provide a value for these. 

```
// INPUTS (OPTIONAL)
```
Workflow inputs which are not required to run the workflow. These are optional. 

```
// PROCESS: ALIGN_AND_TAG
```
Inputs which are specific to a particular process. <br>
These are usually static values rather than files.  <br> 
May be *mandatory* or *optional*. 

<br>

**Setting up Workflow Inputs**

Janis Translate will enter values for workflow inputs where possible. <br>
Others need to be manually supplied as they are specific to the input data you wish to use. 

In our case, we need to supply values for those under the `// INPUTS (MANDATORY)` heading. 
Specifically, we need to provide sample data for the `bams`, `reference`, and `readgroups` inputs. 

Copy and paste the following text, supplying values for these inputs in your current `nextflow.config` file:

```
    // INPUTS (MANDATORY)
    bams        = [
        "../../../../sample_data/cwl/2895499223.bam",
        "../../../../sample_data/cwl/2895499237.bam",
    ]
    reference   = [
        "../../../../sample_data/cwl/chr17_test.fa",
        "../../../../sample_data/cwl/chr17_test.fa.amb",
        "../../../../sample_data/cwl/chr17_test.fa.ann",
        "../../../../sample_data/cwl/chr17_test.fa.bwt",
        "../../../../sample_data/cwl/chr17_test.fa.fai",
        "../../../../sample_data/cwl/chr17_test.dict",
        "../../../../sample_data/cwl/chr17_test.fa.pac",
        "../../../../sample_data/cwl/chr17_test.fa.sa",
    ] 
    readgroups  = [
        '"@RG\tID:2895499223\tPU:H7HY2CCXX.3.ATCACGGT\tSM:H_NJ-HCC1395-HCC1395\tLB:H_NJ-HCC1395-HCC1395-lg24-lib1\tPL:Illumina\tCN:WUGSC"',
        '"@RG\tID:2895499237\tPU:H7HY2CCXX.4.ATCACGGT\tSM:H_NJ-HCC1395-HCC1395\tLB:H_NJ-HCC1395-HCC1395-lg24-lib1\tPL:Illumina\tCN:WUGSC"'
    ]
```

<br>


**Run the Workflow**

Ensure you are in the `translated/` working directory, where `nextflow.config` and `main.nf` reside. 

If not, use the following to change directory. 
```
cd translated/
```

To run the workflow using our sample data, we can now write the following command: 
```
nextflow run main.nf
```

While the workflow runs, you will encounter this error:

```
Access to 'MARK_DUPLICATES_AND_SORT.out' is undefined since the process 'MARK_DUPLICATES_AND_SORT' has not been invoked before accessing the output attribute
```

This is somewhat expected. Janis translate doesn't produce perfect translations - just the best it can do. <br>
This is the first of ***3*** errors we will encounter and fix while making this workflow runnable. 


<br>

## Manual Adjustments

### Error 1: Cutadapt

The first issue we need to address is caused by the cutadapt process translation. 

**Error message**
```
Caused by:
  Missing output file(s) `out1*` expected by process `CUTADAPT (10)`

Command executed:
  cutadapt     MCL1-LD-luminalpregnant.fastq.gz
```

This message tells us that a file called `out1*` should have been present in the working directory, but could not be found. 

>Note<br>
>Nextflow supports ***wildcards*** wherever paths are used. <br>
>Here, `out1*` will collect any file which starts with `out1`, as `*` matches any text. 

<br>

**Diagnosing the Error**

The software manual for a particular tool is a good place to start.
In our case, the [cutadapt user guide](https://cutadapt.readthedocs.io/en/stable/guide.html) tells us what arguments we need to supply. 

From this documentation, we see that the basic usage is as follows
```
cutadapt -a AACCGGTT -o output.fastq input.fastq
```

Further down the cutadapt documentation, we see that the report will be printed to stdout. <br>
We want our script section to resemble the following:
```
cutadapt -a AACCGGTT -o output.fastq input.fastq > report.txt
```

<br>


When diagnosing these types of errors, its always a good idea to look at the process working directory. <br>
Nextflow will provide this path any time a process encounters an error. 

Navigate to this folder. You should see something similar to these files:
```
- .command.begin
- .command.err                      # stderr
- .command.log
- .command.out                      # stdout
- .command.run
- .command.sh                       # the shell command which was executed
- .exitcode
- MCL1-LD-luminalpregnant.fastq.gz
```

We can see that the input `.fastq.gz` file was localised to the working directory, but there is no `out1*` or `report.txt` file which we expect. 

<br>

Open `.command.sh` to view the command which was run. It should be similar to the following:
```
#!/bin/bash -ue
cutadapt     MCL1-LD-luminalpregnant.fastq.gz
```

We see that `cutadapt` is being run on the input file, but we aren't capturing any output. <br>
Open `.command.out` to view stdout. In this case, we see the 'trimmed' reads:
```
@SRR1552447.1 DCV4KXP1:223:C2CTUACXX:1:1101:2099:2233 length=100
GGAAATCCTGAAGAAGATTATTGATATGATCAAGTACATCCAATATCAACAGGTCACCATCCCCCAACTTCCNCAANCTCNTCATCCCCAGATACCTGTG
...
```

<br>

It seems there are a few issues: 
- Our process is missing an input for an adapter sequence to trim
- Our shell command is missing the `-a` and `-o` arguments
- We are not redirecting stdout to a report file


<br>

**Solution**

Update your `CUTADAPT` process definition in `modules/cutadapt.nf`. 

- Add a new `val` process input called `adapters` to accept adapter sequences
- Modify the `script:` section to include the `-a` argument and pass the `${adapters}` input as its value
- Modify the `script:` section to include the `-o` argument and pass `${library_input_1.simpleName}_cutadapt.fastq.gz` as its value
- Modify the `script:` section to pipe stdout to `${library_input_1.simpleName}_cutadapt_report.txt`
- Modify the `outputs:` section so that out1 collects `${library_input_1.simpleName}_cutadapt.fastq.gz`
- Modify the `outputs:` section so that out_report collects `${library_input_1.simpleName}_cutadapt_report.txt`

Your process definition should now look similar to the following: 

```
process CUTADAPT {
    
    container "quay.io/biocontainers/cutadapt:1.16--py35_2"
    publishDir "${params.outdir}/cutadapt"

    input:
    path library_input_1
    val adapter

    output:
    path "${library_input_1.simpleName}_cutadapt.fastq.gz", emit: out1
    path "${library_input_1.simpleName}_cutadapt_report.txt", emit: out_report

    script:
    """
    cutadapt \
    -a ${adapter} \
    -o ${library_input_1.simpleName}_cutadapt.fastq.gz \
    ${library_input_1} \
    > ${library_input_1.simpleName}_cutadapt_report.txt
    """

}
```

<br>

>NOTE <br>
>Why were these `cutadapt` arguments missed? <br><br>
>The default setting for `janis translate` is to translate tools according to their usage in the workflow. <br>
>This is because Galaxy Tool Wrappers are complex and often have many inputs - most of which aren't relevant to the workflow being translated. <br><br>
>To get the extended tool translation we can use `--mode full` when running janis translate:
>```
>janis translate --from galaxy --to nextflow source/rna_seq_reads_to_counts.ga --mode full
>```
>
>For the `CUTADAPT` process, it will now look similar to the following: 
>```
>nextflow.enable.dsl=2
>
>process CUTADAPT {
>    
>    container "quay.io/biocontainers/cutadapt:1.16--py35_2"
>    publishDir "${params.outdir}/cutadapt"
>
>    input:
>    path library_input_1
>    path info_file
>    path output_file
>    path paired_output
>    path rest_file
>    path too_long_output
>    path too_long_paired_output
>    path too_short_output
>    path too_short_paired_output
>    path untrimmed_output
>    path untrimmed_paired_output
>    path wildcard_file
>    val discard
>    val discard_untrimmed
>    val mask_adapter
>    val match_read_wildcards
>    val no_indels
>    val no_trim
>    val trim_n
>    val error_rate
>    val length
>    val length_tag
>    val max_n
>    val maximum_length
>    val minimum_length
>    val nextseq_trim
>    val option_j
>    val option_u
>    val option_u
>    val overlap
>    val pair_filter
>    val prefix
>    val quality_cutoff
>    val read1
>    val read2
>    val suffix
>    val times
>
>    output:
>    path "out1*", emit: out1
>    path "out2*", emit: out2
>    path "info_file.txt", emit: out_info_file
>    path "report.txt", emit: out_report
>    path "rest_output.fastqsanger", emit: out_rest_output
>    path "too_long_output.fastqsanger", emit: out_too_long_output
>    path "too_long_paired_output.fastqsanger", emit: out_too_long_paired_output
>    path "too_short_output.fastqsanger", emit: out_too_short_output
>    path "too_short_paired_output.fastqsanger", emit: out_too_short_paired_output
>    path "untrimmed_output.fastqsanger", emit: out_untrimmed_output
>    path "untrimmed_paired_output.fastqsanger", emit: out_untrimmed_paired_output
>    path "wild_output.txt", emit: out_wild_output
>
>    script:
>    def info_file = info_file.simpleName != params.NULL_VALUE ? "--info-file=${info_file}" : ""
>    def output_file = output_file.simpleName != params.NULL_VALUE ? "--output=${output_file}" : ""
>    def paired_output = paired_output.simpleName != params.NULL_VALUE ? "--paired-output=${paired_output}" : ""
>    def rest_file = rest_file.simpleName != params.NULL_VALUE ? "--rest-file=${rest_file}" : ""
>    def too_long_output = too_long_output.simpleName != params.NULL_VALUE ? "--too-long-output=${too_long_output}" : ""
>    def too_long_paired_output = too_long_paired_output.simpleName != params.NULL_VALUE ? "--too-long-paired-output=${too_long_paired_output}" : ""
>    def too_short_output = too_short_output.simpleName != params.NULL_VALUE ? "--too-short-output=${too_short_output}" : ""
>    def too_short_paired_output = too_short_paired_output.simpleName != params.NULL_VALUE ? "--too-short-paired-output=${too_short_paired_output}" : ""
>    def untrimmed_output = untrimmed_output.simpleName != params.NULL_VALUE ? "--untrimmed-output=${untrimmed_output}" : ""
>    def untrimmed_paired_output = untrimmed_paired_output.simpleName != params.NULL_VALUE ? "--untrimmed-paired-output=${untrimmed_paired_output}" : ""
>    def wildcard_file = wildcard_file.simpleName != params.NULL_VALUE ? "--wildcard-file=${wildcard_file}" : ""
>    def discard = discard ? "--discard" : ""
>    def discard_untrimmed = discard_untrimmed ? "--discard-untrimmed" : ""
>    def mask_adapter = mask_adapter ? "--mask-adapter" : ""
>    def match_read_wildcards = match_read_wildcards ? "--match-read-wildcards" : ""
>    def no_indels = no_indels ? "--no-indels" : ""
>    def no_trim = no_trim ? "--no-trim" : ""
>    def trim_n = trim_n ? "--trim-n" : ""
>    def error_rate = error_rate != params.NULL_VALUE ? error_rate : 0.1
>    def length = length != params.NULL_VALUE ? length : 0
>    def length_tag = length_tag != params.NULL_VALUE ? "--length-tag ${length_tag}" : ""
>    def max_n = max_n != params.NULL_VALUE ? "--max-n ${max_n}" : ""
>    def maximum_length = maximum_length != params.NULL_VALUE ? maximum_length : 0
>    def minimum_length = minimum_length != params.NULL_VALUE ? minimum_length : 0
>    def nextseq_trim = nextseq_trim != params.NULL_VALUE ? nextseq_trim : 0
>    def option_u = option_u != params.NULL_VALUE ? option_u : 0
>    def option_j = option_j != params.NULL_VALUE ? option_j : 1
>    def option_u = option_u != params.NULL_VALUE ? option_u : 0
>    def overlap = overlap != params.NULL_VALUE ? overlap : 3
>    def pair_filter = pair_filter != params.NULL_VALUE ? pair_filter : "any"
>    def prefix = prefix != params.NULL_VALUE ? "--prefix ${prefix}" : ""
>    def quality_cutoff = quality_cutoff != params.NULL_VALUE ? quality_cutoff : "0"
>    def read1 = read1 != params.NULL_VALUE ? read1 : ""
>    def read2 = read2 != params.NULL_VALUE ? read2 : ""
>    def suffix = suffix != params.NULL_VALUE ? "--suffix ${suffix}" : ""
>    def times = times != params.NULL_VALUE ? times : 1
>    """
>    cutadapt \
>    ${library_input_1} \
>    ${info_file} \
>    ${output_file} \
>    ${paired_output} \
>    ${rest_file} \
>    ${too_long_output} \
>    ${too_long_paired_output} \
>    ${too_short_output} \
>    ${too_short_paired_output} \
>    ${untrimmed_output} \
>    ${untrimmed_paired_output} \
>    ${wildcard_file} \
>    --error-rate ${error_rate} \
>    --length ${length} \
>    ${length_tag} \
>    ${max_n} \
>    --maximum-length ${maximum_length} \
>    --minimum-length ${minimum_length} \
>    --nextseq-trim ${nextseq_trim} \
>    --overlap ${overlap} \
>    --pair-filter ${pair_filter} \
>    ${prefix} \
>    --quality-cutoff ${quality_cutoff} \
>    ${suffix} \
>    --times ${times} \
>    -U ${option_u} \
>    -j ${option_j} \
>    -u ${option_u} \
>    ${discard} \
>    ${discard_untrimmed} \
>    ${mask_adapter} \
>    ${match_read_wildcards} \
>    ${no_indels} \
>    ${no_trim} \
>    ${trim_n} \
>    ${read1} \
>    ${read2} \
>    """
>
>}
>```
>
>In this case: 
>- The `--output` option is now present as `--output=${output_file}` when supplied
>- The `--adapter` option is not needed as we are using `--info-file=${info_file}` to supply adapter information


<br>

### Error 2: FastQC Outputs

The second error is due to output collection in the `FASTQC` process. 


**Error message**

```
Caused by:
  Missing output file(s) `output.html` expected by process `FASTQC1 (4)`

Command executed:
  fastqc     MCL1-DJ-basalpregnant.fq.gz
```

This is telling us that `output.html` was expected in the process working directory when collecting outputs, but could not be found. 

<br>

**Diagnosing the Error**

Similar to Error 1, let's have a look at the process working directory to see which files were created. 

It should look similar to the following: 
```
- .
- ..
- .command.begin
- .command.err
- .command.log
- .command.out
- .command.run
- .command.sh
- .exitcode
- MCL1-DJ-basalpregnant.fq.gz -> [path to input data]
- MCL1-DJ-basalpregnant_fastqc.html
- MCL1-DJ-basalpregnant_fastqc.zip
```

Instead of producing `output.html`, we have produced `MCL1-DJ-basalpregnant_fastqc.html`.

It seems that (by default) FastQC produces a `.html` and `.zip` using the input file's base name and adding `_fastqc`.

We will modify the process `outputs:` section to collect these files instead. 

<br>

**Solution**

View the `FASTQC` process in `modules/fastqc.nf`.

We will change the output collection expressions so they capture the correct outputs. 

Modify the `outputs:` section to resemble the following: 
```
    output:
    path "${input_file.simpleName}_fastqc.html", emit: out_html_file
    path "${input_file.simpleName}_fastqc.zip", emit: out_text_file
```

<br>

### Error 3: Hisat2 Inputs

The 3rd error is due to difficulties in translating some Galaxy Tool Wrappers. 


**Error message**

```
Caused by:
  Process `HISAT2 (1)` terminated with an error exit status (1)

Command executed:
  hisat2     MCL1-DK-basallactate_cutadapt.fastq.gz
```


<br>

**Diagnosing the Error**

As in the previous errors, let's have a look at the process working directory. 

This time, we will look at `.command.err` which is the ***stderr*** of our command. Most command line software will print error messages to ***stderr***. 

The first lines of `.command.err` are hisat2-align telling us the command line format to run this software: 

```
Usage: 
  hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
```
At the bottom of `.command.err` we see hisat2-align telling us what we missed in our command: 
```
Error: Must specify at least one read input with -U/-1/-2
(ERR): hisat2-align exited with value 1
```

From this information, we can gather that the `script:` section of our `HISAT2` process isn't formatted correctly. 

<br>

View the `HISAT2` process definition in `modules/hisat2.nf`. 
It should look similar to the following: 
```
process HISAT2 {
    
    container "quay.io/biocontainers/hisat2:2.2.1--h87f3376_5"
    publishDir "${params.outdir}/hisat2"

    input:
    path library_input_1

    output:
    path "summary.txt", emit: out_summary_file
    path "unknown_collection_pattern", emit: output_alignments

    script:
    """
    hisat2 \
    ${library_input_1} \
    """

}
```

The process definition has a few errors. 
- We need a new process input for the hisat2 index we will align against.
- We need to supply the hisat2 index using `-x` in the script section.
- As our reads are single-end, we need to supply these using `-U`in the script section. 

The is also an error with the `outputs:` section of the `HISAT2` process, but we will tackle this in the next error. 

>Note<br>
>Why didn't `janis translate` create a process input for the hisat2 index?<br><br>
>Galaxy Tool Wrappers can use inbuilt data on the Galaxy server when running tools. <br>
>In this case, the `hisat2` Galaxy Tool Wrapper uses prebuilt indexes which you can select using a dropdown.<br>
>As this happens internally and isn't part of the canonical `hisat2` software command, `janis translate` doesn't translate it. <br><br>
>In the future, `janis translate` may support inbuilt data. <br>
>Care needs to be taken due to the sheer size of these files.

<br>

**Solution**

modules/hisat2.nf
```
process HISAT2 {
    
    container "quay.io/biocontainers/hisat2:2.2.1--h87f3376_5"
    publishDir "${params.outdir}/hisat2"

    input:
    path library_input_1
    path index

    output:
    path "summary.txt", emit: out_summary_file
    path "unknown_collection_pattern", emit: output_alignments

    script:
    """
    hisat2 \
    -x ${index[0].simpleName} \
    -U ${library_input_1} \
    """

}
```

nextflow.config
```
    ...
    multiqc_config  = "templates/multiqc_config" 
    adapter         = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    hisat2_index    = "(local_dir)/sample_data/galaxy/rnaseq_reads_to_counts_workflow/hisat2_index/*.ht2"
```

main.nf
```
// data which will be passed as channels
ch_collection_column_join_script  = Channel.fromPath( params.collection_column_join_script )
ch_in_input_fastqs_collection     = Channel.fromPath( params.in_input_fastqs_collection ).toList()
ch_in_input_reference_gene_bed    = Channel.fromPath( params.in_input_reference_gene_bed )
ch_multiqc_config                 = Channel.fromPath( params.multiqc_config )
ch_hisat2_index                   = Channel.fromPath( params.hisat2_index ).toList()
```


### Error 4: Hisat2 Outputs



