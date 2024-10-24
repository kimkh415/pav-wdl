version 1.0

task call_inv_flag_insdel_cluster_indel_hap {
  input {
    String sample
    File pav_conf
    File pav_asm
    File inbed
    String hap
    String threads
    String mem_gb
    String vartype
  }
  String docker_dir = "/assemblybased"
  String work_dir = "/cromwell_root/assemblybased"
  command <<<
    set -euxo pipefail
    cd ~{docker_dir}
    source activate lr-pav
    mkdir -p ~{work_dir}
    cd ~{work_dir}
    cp -r ~{docker_dir}/pav .
    
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{inbed}
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/inv_caller/flag/insdel_~{vartype}_~{hap}.bed.gz
    tar zcvf call_inv_flag_insdel_cluster_indel_~{sample}_~{vartype}_~{hap}.tgz temp/~{sample}/inv_caller/flag/insdel_~{vartype}_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(work_dir + "/.snakemake/log/*.snakemake.log")
    File bed = work_dir + "/call_inv_flag_insdel_cluster_indel_~{sample}_~{vartype}_~{hap}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "fcunial/assemblybased"
  }
}

task call_inv_flag_insdel_cluster_sv_hap {
  input {
    String sample
    File pav_conf
    File pav_asm
    File inbed
    String hap
    String threads
    String mem_gb
    String vartype
  }
  String docker_dir = "/assemblybased"
  String work_dir = "/cromwell_root/assemblybased"
  command <<<
    set -euxo pipefail
    cd ~{docker_dir}
    source activate lr-pav
    mkdir -p ~{work_dir}
    cd ~{work_dir}
    cp -r ~{docker_dir}/pav .
    
    cp ~{pav_conf} ./config.json
    tar zxvf ~{inbed}
    tar zxvf ~{pav_asm}
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/inv_caller/flag/insdel_~{vartype}_~{hap}.bed.gz
    tar zcvf call_inv_flag_insdel_cluster_indel_~{sample}_~{vartype}_~{hap}.tgz temp/~{sample}/inv_caller/flag/insdel_~{vartype}_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(work_dir + "/.snakemake/log/*.snakemake.log")
    File bed = work_dir + "/call_inv_flag_insdel_cluster_indel_~{sample}_~{vartype}_~{hap}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "fcunial/assemblybased"
  }
}

task call_inv_cluster_indel_hap {
  input {
    String sample
    File pav_conf
    File pav_asm
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  }
  String docker_dir = "/assemblybased"
  String work_dir = "/cromwell_root/assemblybased"
  command <<<
    set -euxo pipefail
    cd ~{docker_dir}
    source activate lr-pav
    mkdir -p ~{work_dir}
    cd ~{work_dir}
    cp -r ~{docker_dir}/pav .
    
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{inbed}
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/inv_caller/flag/cluster_~{vartype}_~{hap}.bed.gz
    tar zcvf call_inv_cluster_indel_~{hap}_~{sample}_~{vartype}.tgz temp/~{sample}/inv_caller/flag/cluster_~{vartype}_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(work_dir + "/.snakemake/log/*.snakemake.log")
    File bed = work_dir + "/call_inv_cluster_indel_~{hap}_~{sample}_~{vartype}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "fcunial/assemblybased"
  }
}

task call_inv_cluster_snv_hap {
  input {
    String sample
    File pav_conf
    File pav_asm
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  }
  String docker_dir = "/assemblybased"
  String work_dir = "/cromwell_root/assemblybased"
  command <<<
    set -euxo pipefail
    cd ~{docker_dir}
    source activate lr-pav
    mkdir -p ~{work_dir}
    cd ~{work_dir}
    cp -r ~{docker_dir}/pav .
  
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{inbed}
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/inv_caller/flag/cluster_~{vartype}_~{hap}.bed.gz
    tar zcvf call_inv_cluster_snv_~{hap}_~{sample}_~{vartype}.tgz temp/~{sample}/inv_caller/flag/cluster_~{vartype}_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(work_dir + "/.snakemake/log/*.snakemake.log")
    File bed = work_dir + "/call_inv_cluster_snv_~{hap}_~{sample}_~{vartype}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "fcunial/assemblybased"
  }
}

task call_inv_merge_flagged_loci_hap {
  input {
    String sample
    String hap
    File pav_conf
    File pav_asm
    File indelFlag
    File svFlag
    File snvCluster
    File indelCluster
    String threads
    String mem_gb
  }
  String docker_dir = "/assemblybased"
  String work_dir = "/cromwell_root/assemblybased"
  command <<<
    set -euxo pipefail
    cd ~{docker_dir}
    source activate lr-pav
    mkdir -p ~{work_dir}
    cd ~{work_dir}
    cp -r ~{docker_dir}/pav .
    
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{indelFlag}
    tar zxvf ~{svFlag}
    tar zxvf ~{snvCluster}
    tar zxvf ~{indelCluster}
    tree
    snakemake -s pav/Snakefile --cores ~{threads} results/~{sample}/inv_caller/flagged_regions_~{hap}.bed.gz
    tar zcvf call_inv_merge_flagged_loci_~{hap}_~{sample}.tgz results/~{sample}/inv_caller/flagged_regions_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(work_dir + "/.snakemake/log/*.snakemake.log")
    File bed = work_dir + "/call_inv_merge_flagged_loci_~{hap}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "fcunial/assemblybased"
  }
}

task call_inv_batch_hap {
  input {
    String sample
    File pav_conf
    File pav_asm
    String hap
    File trimBed
    File flag
    File asmGz
    File refGz
    String batch
    String threads
    String mem_gb
  }
  String docker_dir = "/assemblybased"
  String work_dir = "/cromwell_root/assemblybased"
  command <<<
    set -euxo pipefail
    cd ~{docker_dir}
    source activate lr-pav
    mkdir -p ~{work_dir}
    cd ~{work_dir}
    cp -r ~{docker_dir}/pav .

    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{trimBed}
    tar zxvf ~{flag}
    tar zxvf ~{asmGz}
    tar zxvf ~{refGz}
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/inv_caller/batch/~{hap}/inv_call_~{batch}.bed.gz
    tar zcvf call_inv_batch_~{hap}_~{batch}_~{sample}.tgz temp/~{sample}/inv_caller/batch/~{hap}/inv_call_~{batch}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(work_dir + "/.snakemake/log/*.snakemake.log")
    File bed = work_dir + "/call_inv_batch_~{hap}_~{batch}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "fcunial/assemblybased"
  }
}

task call_inv_batch_merge_hap {
  input {
    String sample
    File pav_conf
    File pav_asm
    String hap
    Array[File] invBed
    String threads
    String mem_gb
  }
  String docker_dir = "/assemblybased"
  String work_dir = "/cromwell_root/assemblybased"
  command <<<
    set -euxo pipefail
    cd ~{docker_dir}
    source activate lr-pav
    mkdir -p ~{work_dir}
    cd ~{work_dir}
    cp -r ~{docker_dir}/pav .

    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    echo ~{sep=" " invBed} | tr " " "\n" | xargs -I '@' tar zxvf @
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/inv_caller/sv_inv_~{hap}.bed.gz
    tar zcvf call_inv_batch_merge_~{hap}_~{sample}.tgz temp/~{sample}/inv_caller/sv_inv_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(work_dir + "/.snakemake/log/*.snakemake.log")
    File bed = work_dir + "/call_inv_batch_merge_~{hap}_~{sample}.tgz "
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "fcunial/assemblybased"
  }
}