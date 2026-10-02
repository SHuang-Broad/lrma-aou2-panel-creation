version 1.0

workflow SubsetAndSplitVcf {

    input {
        File joint_vcf
        File joint_vcf_idx
        String region
        String gcs_output_dir
        String output_tag           # per-sample BCFs will be copied to gcs_output_dir/{sample_name}.{output_tag}.bcf

        # optional hierarchical split
        File? sample_batches_tsv      # see bcftools +split --help; e.g., contains rows: {sample_name_1},{sample_name_2},...\t-\tbatch-0
    }

    call SubsetAndSplitVcf { input:
        vcf = joint_vcf,
        vcf_idx = joint_vcf_idx,
        region = region,
        gcs_output_dir = if defined(sample_batches_tsv) then gcs_output_dir + "/batches" else gcs_output_dir,
        sample_batches_tsv = sample_batches_tsv,
        output_tag = output_tag
    }
    
    if (defined(sample_batches_tsv)) {
        Int num_batches = length(read_lines(select_first([sample_batches_tsv])))
        scatter (i in range(num_batches)) {
            call SubsetAndSplitVcf as SubsetAndSplitVcfBatch { input:
                vcf = gcs_output_dir + "/batches/batch-" + i + "." + output_tag + ".bcf",
                gcs_output_dir = gcs_output_dir,
                wait_for_me = SubsetAndSplitVcf.number_of_sites,
                output_tag = output_tag
            }
        }
    }

    output {
        Int number_of_sites = SubsetAndSplitVcf.number_of_sites
        Array[String] split_vcf_paths = select_first([flatten(select_first([SubsetAndSplitVcfBatch.split_vcf_paths])), 
                                                      SubsetAndSplitVcf.split_vcf_paths])
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Boolean? use_ssd
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

task SubsetAndSplitVcf {

    input {
        File vcf
        File? vcf_idx       # only needed for initial stream, not for batches
        String? region      # only needed for initial stream, not for batches
        String gcs_output_dir
        String output_tag
        File? sample_batches_tsv
        Int? wait_for_me    # dummy input to gate scatter until initial split is done
        Int view_verbosity = 8
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        vcf: {
            localization_optional: true
        }
    }
    
    String view_input_arg = if defined(vcf_idx) then "\"~{vcf}##idx##~{vcf_idx}\"" else "~{vcf}"

    command <<<
        set -euxo pipefail

        # see https://github.com/samtools/htslib/issues/803#issuecomment-444514336, https://github.com/broadinstitute/bcftools-patched
        mkfifo /tmp/token_fifo
        ( while true ; do curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        export HTS_AUTH_LOCATION="/tmp/token_fifo"

        # we stream to an intermediate file, since piping directly to bcftools +split still results in GOAWAY/Libcurl issues
        bcftools view --no-version ~{view_input_arg} \
            ~{"--regions " + region} \
            --regions-overlap 0 \
            --verbosity ~{view_verbosity} \
            -Ob -o ~{output_tag}.bcf

        # check number of sites to guard against streaming errors
        bcftools view -H ~{output_tag}.bcf | wc -l > number_of_sites.txt

        bcftools +split ~{output_tag}.bcf ~{"--samples-file " + sample_batches_tsv} -Ob -o output

        for bcf in output/*.bcf; do
            bcf_basename=$(basename $bcf .bcf)
            mv $bcf output/$bcf_basename.~{output_tag}.bcf
        done

        gcloud storage cp output/*.bcf ~{gcs_output_dir}/
        gsutil ls ~{gcs_output_dir}/*bcf > output_vcf_paths.txt
    >>>

    output {
        Array[String] split_vcf_paths = read_lines("output_vcf_paths.txt")
        Int number_of_sites = read_int("number_of_sites.txt")
    }
    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            20,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  5,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsde-methods/slee/lr-bcftools-patched-gcloud/lr-bcftools-patched-gcloud:1.23"      # see https://github.com/broadinstitute/bcftools-patched
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
