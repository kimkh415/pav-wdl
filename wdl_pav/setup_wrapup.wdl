version 1.0

task tar_asm {
  input {
    File ref
    File hapOne
    File hapTwo
    String sample
    String threads
    String mem_gb
  }
  command <<<
    set -eux
    mkdir -p asm/~{sample}
    cp ~{ref} asm/ref.fa
    samtools faidx asm/ref.fa
    cp ~{hapOne} asm/~{sample}/h1.fa.gz
    cp ~{hapTwo} asm/~{sample}/h2.fa.gz
    tar zcvf asm.tgz asm/
  >>>
  output {
    File asm_tar = "asm.tgz"
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

task call_final_bed {
  input {
    File pav_conf
    File pav_asm
    File invBed
    File insBed
    File delBed
    File snvBed
    String threads
    String mem_gb
    String sample
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
    tar zxvf ~{invBed}
    tar zxvf ~{snvBed}
    tar zxvf ~{insBed}
    tar zxvf ~{delBed}
    tree
    snakemake -s pav/Snakefile --cores ~{threads} results/~{sample}/bed/snv_snv.bed.gz results/~{sample}/bed/indel_ins.bed.gz results/~{sample}/bed/indel_del.bed.gz results/~{sample}/bed/sv_ins.bed.gz results/~{sample}/bed/sv_del.bed.gz results/~{sample}/bed/sv_inv.bed.gz results/~{sample}/bed/fa/indel_ins.fa.gz results/~{sample}/bed/fa/indel_del.fa.gz results/~{sample}/bed/fa/sv_ins.fa.gz results/~{sample}/bed/fa/sv_del.fa.gz results/~{sample}/bed/fa/sv_inv.fa.gz
    tar zcvf final_bed.tgz results/~{sample}/bed/snv_snv.bed.gz results/~{sample}/bed/indel_ins.bed.gz results/~{sample}/bed/indel_del.bed.gz results/~{sample}/bed/sv_ins.bed.gz results/~{sample}/bed/sv_del.bed.gz results/~{sample}/bed/sv_inv.bed.gz results/~{sample}/bed/fa/indel_ins.fa.gz results/~{sample}/bed/fa/indel_del.fa.gz results/~{sample}/bed/fa/sv_ins.fa.gz results/~{sample}/bed/fa/sv_del.fa.gz results/~{sample}/bed/fa/sv_inv.fa.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(work_dir + "/.snakemake/log/*.snakemake.log")
    File bed = work_dir + "/final_bed.tgz"
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

task data_ref_contig_table{
  input {
    File pav_conf
    File pav_asm
    String threads
    String mem_gb
    String sample
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
    tree
    snakemake -s snakemake -s pav/Snakefile --cores ~{threads} data/ref/contig_info.tsv.gz
    tar zcvf contig_info.tgz data/ref/contig_info.tsv.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(work_dir + "/.snakemake/log/*.snakemake.log")
    File contigInfo = work_dir + "/contig_info.tgz"
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

task write_vcf {
  input {
    File pav_conf
    File pav_asm
    File contigInfo
    File finalBedOut
    String threads
    String mem_gb
    String sample
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
    tar zxvf ~{contigInfo}
    tar zxvf ~{finalBedOut}
    tree
    snakemake -s pav/Snakefile --cores ~{threads} pav_~{sample}.vcf.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(work_dir + "/.snakemake/log/*.snakemake.log")
    File vcf = work_dir + "/pav_~{sample}.vcf.gz"
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