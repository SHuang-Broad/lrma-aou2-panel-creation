version 1.0

workflow SubsetAndSplitVcf {

    input {
        File joint_vcf
        File joint_vcf_idx
        String region
        Int expected_number_of_sites

        String gcs_output_dir
        String output_tag

        # optional hierarchical split
        File? sample_batches_tsv
    }

    output {
        Int number_of_sites = SubsetAndSplitVcf.number_of_sites
        Array[String] split_vcf_paths = select_first([flatten(select_first([SubsetAndSplitVcfBatch.split_vcf_paths])),
                                                      SubsetAndSplitVcf.split_vcf_paths])
    }

    parameter_meta {
        output_tag: {
            description: "per-sample BCFs will be copied to gcs_output_dir/{sample_name}.{output_tag}.bcf"
        }
        sample_batches_tsv: {
            description: "see bcftools +split --help; e.g., contains rows: {sample_name_1},{sample_name_2},...\t-\tbatch-0"
        }
        expected_number_of_sites: {
            description: "expected number of sites in this region, used for safeguard against streaming errors"
        }
    }

    call SubsetAndSplitVcf { input:
        vcf     = joint_vcf,
        vcf_idx = joint_vcf_idx,

        region  = region,
        expected_number_of_sites = expected_number_of_sites,

        sample_batches_tsv = sample_batches_tsv,

        gcs_output_dir = if defined(sample_batches_tsv) then gcs_output_dir + "/batches" else gcs_output_dir,
        output_tag     = output_tag
    }
    
    if (defined(sample_batches_tsv)) {
        scatter (batch_vcf in SubsetAndSplitVcf.split_vcf_paths) {
            call SubsetAndSplitVcf as SubsetAndSplitVcfBatch { input:
                vcf = batch_vcf,

                expected_number_of_sites = expected_number_of_sites,

                gcs_output_dir = gcs_output_dir,
                output_tag = output_tag
            }
        }
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
        File  vcf
        File? vcf_idx       # only needed for initial stream, not for batches

        Int expected_number_of_sites

        String? region      # only needed for initial stream, not for batches

        File? sample_batches_tsv

        String gcs_output_dir
        String output_tag

        Int view_verbosity = 8

        RuntimeAttr? runtime_attr_override
    }

    output {
        Int number_of_sites = read_int("number_of_sites.txt")
        Array[String] split_vcf_paths = read_lines("output_vcf_paths.txt")
    }

    parameter_meta {
        vcf: {
            localization_optional: true
        }

        sample_batches_tsv: {
            description: ""
        }

        split_vcf_paths: {
            description: ""
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
        N=$(bcftools view -H ~{output_tag}.bcf | wc -l | awk '{print $1}')
        if [[ ${N} -ne ~{expected_number_of_sites} ]]; then echo "Number of sites unexpected"; exit 1; fi

        bcftools +split ~{output_tag}.bcf ~{"--samples-file " + sample_batches_tsv} -Ob -o output

        for bcf in output/*.bcf; do
            bcf_basename=$(basename $bcf .bcf)
            mv $bcf output/$bcf_basename.~{output_tag}.bcf
        done

        # user-controlled result saving
        gcloud storage cp \
            output/*.bcf \
            "~{gcs_output_dir}/"

        # Report only the files produced by this invocation.
        # Do not list the shared GCS directory: other scattered calls may be
        # writing there concurrently, and previous runs may have left files.
        : > output_vcf_paths.txt
        for bcf in output/*.bcf; do
            bcf_basename=$(basename "$bcf")
            echo "~{gcs_output_dir}/$bcf_basename" >> output_vcf_paths.txt
        done
    >>>

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
