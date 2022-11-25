version 1.0

import "wdl_pav/setup_wrapup.wdl" as setup
import "wdl_pav/align.wdl" as align
import "wdl_pav/call.wdl" as call_pav
import "wdl_pav/call_inv.wdl" as call_inv
import "wdl_pav/call_lg.wdl" as call_lg


workflow pav {
  input {
    File ref
    File refFai
    File hapOne
    File hapTwo
    String sample
    File config
  }
  parameter_meta {
      ref: "In .fa format."
      hapOne: "In .fa.gz format."
      hapTwo: "In .fa.gz format."
      config: "Default, fixed config file provided with this workflow. Contains no information about the sample."
  }

  Array[Array[String]] chroms = read_tsv(refFai)

  call setup.tar_asm {
    input:
      ref = ref,
      hapOne = hapOne,
      hapTwo = hapTwo,
      threads = "1",
      mem_gb = "2",
      sample = sample
  }
  call align.align_ref {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call align.align_get_tig_fa_hap as align_get_tig_fa_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call align.align_get_tig_fa_hap as align_get_tig_fa_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call align.align_ref_anno_n_gap {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      ref_gz = align_ref.refGz,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call align.align_map_hap as align_map_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      refGz = align_ref.refGz,
      asmGz = align_get_tig_fa_h1.asmGz,
      threads = "8",
      mem_gb = "48",
      sample = sample
  }
  call align.align_map_hap as align_map_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      refGz = align_ref.refGz,
      asmGz = align_get_tig_fa_h2.asmGz,
      threads = "8",
      mem_gb = "48",
      sample = sample
  }
  call align.align_get_read_bed_hap as align_get_read_bed_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      refGz = align_ref.refGz,
      tigFa = align_get_tig_fa_h1.asmGz,
      samGz = align_map_h1.samGz,
      threads = "6",
      mem_gb = "32",
      sample = sample
  }
  call align.align_get_read_bed_hap as align_get_read_bed_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      refGz = align_ref.refGz,
      hap = "h2",
      tigFa = align_get_tig_fa_h2.asmGz,
      samGz = align_map_h2.samGz,
      threads = "6",
      mem_gb = "32",
      sample = sample
  }
  call align.align_cut_tig_overlap_hap as align_cut_tig_overlap_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      bedGz = align_get_read_bed_h1.bedGz,
      asmGz = align_get_tig_fa_h1.asmGz,
      threads = "2",
      mem_gb = "8",
      sample = sample
  }
  call align.align_cut_tig_overlap_hap as align_cut_tig_overlap_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      bedGz = align_get_read_bed_h2.bedGz,
      asmGz = align_get_tig_fa_h2.asmGz,
      threads = "2",
      mem_gb = "8",
      sample = sample
  }
  scatter(i in range(10)) {
     call call_pav.call_cigar_hap as call_cigar_h1 {
      input:
        pav_conf = config,
        pav_asm = tar_asm.asm_tar,
        refGz = align_ref.refGz,
        hap = "h1",
        trimBed = align_cut_tig_overlap_h1.trimBed,
        asmGz = align_get_tig_fa_h1.asmGz,
        batch = i,
        threads = "1",
        mem_gb = "24",
        sample = sample
     }
     call call_pav.call_cigar_hap as call_cigar_h2 {
      input:
        pav_conf = config,
        refGz = align_ref.refGz,
        pav_asm = tar_asm.asm_tar,
        hap = "h2",
        trimBed = align_cut_tig_overlap_h2.trimBed,
        asmGz = align_get_tig_fa_h2.asmGz,
        batch = i,
        threads = "1",
        mem_gb = "24",
        sample = sample
     }
  }
  call call_lg.call_lg_split_hap as call_lg_split_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      trimBed = align_cut_tig_overlap_h1.trimBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_lg.call_lg_split_hap as call_lg_split_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      trimBed = align_cut_tig_overlap_h2.trimBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  scatter(i in range(10)) {
     call call_lg.call_lg_discover_hap as call_lg_discover_h1 {
      input:
        pav_conf = config,
        pav_asm = tar_asm.asm_tar,
        refGz = align_ref.refGz,
        hap = "h1",
        trimBed = align_cut_tig_overlap_h1.trimBed,
        gaps = align_ref_anno_n_gap.gaps,
        asmGz = align_get_tig_fa_h1.asmGz,
        batchFile = call_lg_split_h1.batch,
        batch = i,
        threads = "8",
        mem_gb = "32",
        sample = sample
     }
     call call_lg.call_lg_discover_hap as call_lg_discover_h2 {
      input:
        pav_conf = config,
        pav_asm = tar_asm.asm_tar,
        refGz = align_ref.refGz,
        hap = "h2",
        trimBed = align_cut_tig_overlap_h2.trimBed,
        gaps = align_ref_anno_n_gap.gaps,
        asmGz = align_get_tig_fa_h2.asmGz,
        batchFile = call_lg_split_h2.batch,
        batch = i,
        threads = "8",
        mem_gb = "32",
        sample = sample
     }
  }
  call call_pav.call_cigar_merge_hap as call_cigar_merge_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      snvBatch = call_cigar_h1.snvBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_cigar_merge_hap as call_cigar_merge_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      snvBatch = call_cigar_h2.snvBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_lg.call_merge_lg_del_hap as call_merge_lg_del_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      svtype = "del",
      inbed = call_lg_discover_h1.allBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_lg.call_merge_lg_ins_hap as call_merge_lg_ins_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      svtype = "ins",
      inbed = call_lg_discover_h1.allBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_lg.call_merge_lg_inv_hap as call_merge_lg_inv_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      svtype = "inv",
      inbed = call_lg_discover_h1.allBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_lg.call_merge_lg_del_hap as call_merge_lg_del_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      svtype = "del",
      inbed = call_lg_discover_h2.allBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_lg.call_merge_lg_ins_hap as call_merge_lg_ins_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      svtype = "ins",
      inbed = call_lg_discover_h2.allBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_lg.call_merge_lg_inv_hap as call_merge_lg_inv_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      svtype = "inv",
      inbed = call_lg_discover_h2.allBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_inv.call_inv_cluster_indel_hap as call_inv_cluster_indel_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      vartype="indel",
      inbed = call_cigar_merge_h1.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_cluster_snv_hap as call_inv_cluster_snv_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      vartype="snv",
      inbed = call_cigar_merge_h1.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_cluster_indel_hap as call_inv_cluster_indel_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      vartype="indel",
      inbed = call_cigar_merge_h2.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_cluster_snv_hap as call_inv_cluster_snv_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      vartype="snv",
      inbed = call_cigar_merge_h2.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_flag_insdel_cluster_indel_hap as call_inv_flag_insdel_cluster_indel_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      vartype="indel",
      inbed = call_cigar_merge_h1.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_flag_insdel_cluster_indel_hap as call_inv_flag_insdel_cluster_indel_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      vartype="indel",
      inbed = call_cigar_merge_h2.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_flag_insdel_cluster_sv_hap as call_inv_flag_insdel_cluster_sv_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      vartype="sv",
      inbed = call_cigar_merge_h1.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_flag_insdel_cluster_sv_hap as call_inv_flag_insdel_cluster_sv_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      vartype="sv",
      inbed = call_cigar_merge_h2.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_pav.call_mappable_bed_hap as call_mappable_bed_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      trimBed = align_cut_tig_overlap_h1.trimBed,
      delBed = call_merge_lg_del_h1.mergeBed,
      insBed = call_merge_lg_ins_h1.mergeBed,
      invBed = call_merge_lg_inv_h1.mergeBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_mappable_bed_hap as call_mappable_bed_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      trimBed = align_cut_tig_overlap_h2.trimBed,
      delBed = call_merge_lg_del_h2.mergeBed,
      insBed = call_merge_lg_ins_h2.mergeBed,
      invBed = call_merge_lg_inv_h2.mergeBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_inv.call_inv_merge_flagged_loci_hap as call_inv_merge_flagged_loci_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      indelFlag = call_inv_flag_insdel_cluster_indel_h1.bed,
      svFlag = call_inv_flag_insdel_cluster_sv_h1.bed,
      snvCluster = call_inv_cluster_snv_h1.bed,
      indelCluster = call_inv_cluster_indel_h1.bed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_inv.call_inv_merge_flagged_loci_hap as call_inv_merge_flagged_loci_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      indelFlag = call_inv_flag_insdel_cluster_indel_h2.bed,
      svFlag = call_inv_flag_insdel_cluster_sv_h2.bed,
      snvCluster = call_inv_cluster_snv_h2.bed,
      indelCluster = call_inv_cluster_indel_h2.bed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  scatter(i in range(60)) {
     call call_inv.call_inv_batch_hap as call_inv_batch_h1 {
      input:
        pav_conf = config,
        pav_asm = tar_asm.asm_tar,
        hap = "h1",
        trimBed = align_cut_tig_overlap_h1.trimBed,
        flag = call_inv_merge_flagged_loci_h1.bed,
        asmGz = align_get_tig_fa_h1.asmGz,
        batch = i,
        refGz = align_ref.refGz,
        threads = "8",
        mem_gb = "8",
        sample = sample
     }
     call call_inv.call_inv_batch_hap as call_inv_batch_h2 {
      input:
        pav_conf = config,
        pav_asm = tar_asm.asm_tar,
        hap = "h2",
        trimBed = align_cut_tig_overlap_h2.trimBed,
        flag = call_inv_merge_flagged_loci_h2.bed,
        asmGz = align_get_tig_fa_h2.asmGz,
        batch = i,
        refGz = align_ref.refGz,
        threads = "8",
        mem_gb = "8",
        sample = sample
     }
  }
  call call_inv.call_inv_batch_merge_hap as call_inv_batch_merge_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      invBed = call_inv_batch_h1.bed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_inv.call_inv_batch_merge_hap as call_inv_batch_merge_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      invBed = call_inv_batch_h2.bed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_integrate_sources_hap as call_integrate_sources_h1 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h1",
      preInvSvBed = call_cigar_merge_h1.insdelBedMerge,
      delBedIn = call_merge_lg_del_h1.mergeBed,
      insBedIn = call_merge_lg_ins_h1.mergeBed,
      invBedIn = call_merge_lg_inv_h1.mergeBed,
      invBatch = call_inv_batch_merge_h1.bed,
      threads = "2",
      mem_gb = "32",
      sample = sample
  }
  call call_pav.call_integrate_sources_hap as call_integrate_sources_h2 {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      hap = "h2",
      preInvSvBed = call_cigar_merge_h2.insdelBedMerge,
      delBedIn = call_merge_lg_del_h2.mergeBed,
      insBedIn = call_merge_lg_ins_h2.mergeBed,
      invBedIn = call_merge_lg_inv_h2.mergeBed,
      invBatch = call_inv_batch_merge_h2.bed,
      threads = "2",
      mem_gb = "32",
      sample = sample
  }
  scatter(chrom in chroms) {
     call call_pav.call_merge_haplotypes_chrom_svindel as call_merge_haplotypes_chrom_svindel_ins {
      input:
        pav_conf = config,
        pav_asm = tar_asm.asm_tar,
        svtype = "svindel_ins",
        svindel_bed_h1 = call_integrate_sources_h1.all_vars_bed,
        svindel_bed_h2 = call_integrate_sources_h2.all_vars_bed,
        callable_h1 = call_mappable_bed_h2.bed,
        callable_h2 = call_mappable_bed_h1.bed,
        chrom = chrom[0],
        threads = "8",
        mem_gb = "12",
        sample = sample
     }
     call call_pav.call_merge_haplotypes_chrom_svindel as call_merge_haplotypes_chrom_svindel_del {
      input:
        pav_conf = config,
        pav_asm = tar_asm.asm_tar,
        svtype = "svindel_del",
        svindel_bed_h1 = call_integrate_sources_h1.all_vars_bed,
        svindel_bed_h2 = call_integrate_sources_h2.all_vars_bed,
        callable_h1 = call_mappable_bed_h2.bed,
        callable_h2 = call_mappable_bed_h1.bed,
        chrom = chrom[0],
        threads = "8",
        mem_gb = "12",
        sample = sample
     }
     call call_pav.call_merge_haplotypes_chrom_svindel as call_merge_haplotypes_chrom_svinv {
      input:
        pav_conf = config,
        pav_asm = tar_asm.asm_tar,
        svtype = "sv_inv",
        svindel_bed_h1 = call_integrate_sources_h1.all_vars_bed,
        svindel_bed_h2 = call_integrate_sources_h2.all_vars_bed,
        callable_h1 = call_mappable_bed_h2.bed,
        callable_h2 = call_mappable_bed_h1.bed,
        chrom = chrom[0],
        threads = "8",
        mem_gb = "12",
        sample = sample
     }
  }
  scatter(chrom in chroms) {
     call call_pav.call_merge_haplotypes_chrom_snv {
      input:
        pav_conf = config,
        pav_asm = tar_asm.asm_tar,
        svtype = "snv_snv",
        snvBed_h1 = call_integrate_sources_h1.all_vars_bed,
        snvBed_h2 = call_integrate_sources_h2.all_vars_bed,
        callable_h1 = call_mappable_bed_h2.bed,
        callable_h2 = call_mappable_bed_h1.bed,
        chrom = chrom[0],
        threads = "8",
        mem_gb = "12",
        sample = sample
     }
  }
  call call_pav.call_merge_haplotypes as call_merge_haplotypes_snv {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      svtype = "snv_snv",
      inbed = call_merge_haplotypes_chrom_snv.bed,
      callable_h1 = call_mappable_bed_h2.bed,
      callable_h2 = call_mappable_bed_h1.bed,
      integrated_h1 = call_integrate_sources_h1.all_vars_bed,
      integrated_h2 = call_integrate_sources_h2.all_vars_bed,
      threads = "12",
      mem_gb = "24",
      sample = sample
  }
  call call_pav.call_merge_haplotypes as call_merge_haplotypes_inv {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      svtype = "sv_inv",
      inbed = call_merge_haplotypes_chrom_svinv.bed,
      callable_h1 = call_mappable_bed_h2.bed,
      callable_h2 = call_mappable_bed_h1.bed,
      integrated_h1 = call_integrate_sources_h1.all_vars_bed,
      integrated_h2 = call_integrate_sources_h2.all_vars_bed,
      threads = "12",
      mem_gb = "24",
      sample = sample
  }
  call call_pav.call_merge_haplotypes as call_merge_haplotypes_svindel_ins {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      integrated_h1 = call_integrate_sources_h1.all_vars_bed,
      integrated_h2 = call_integrate_sources_h2.all_vars_bed,
      callable_h1 = call_mappable_bed_h2.bed,
      callable_h2 = call_mappable_bed_h1.bed,
      svtype = "svindel_ins",
      inbed = call_merge_haplotypes_chrom_svindel_ins.bed,
      threads = "12",
      mem_gb = "24",
      sample = sample
  }
  call call_pav.call_merge_haplotypes as call_merge_haplotypes_svindel_del {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      svtype = "svindel_del",
      callable_h1 = call_mappable_bed_h2.bed,
      callable_h2 = call_mappable_bed_h1.bed,
      inbed = call_merge_haplotypes_chrom_svindel_del.bed,
      integrated_h1 = call_integrate_sources_h1.all_vars_bed,
      integrated_h2 = call_integrate_sources_h2.all_vars_bed,
      threads = "12",
      mem_gb = "24",
      sample = sample
  }
  call setup.call_final_bed {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      invBed = call_merge_haplotypes_inv.bed,
      insBed = call_merge_haplotypes_svindel_ins.bed,
      delBed = call_merge_haplotypes_svindel_del.bed,
      snvBed = call_merge_haplotypes_snv.bed,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call setup.data_ref_contig_table {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      threads = "1",
      mem_gb = "4",
      sample = sample
  }
  call setup.write_vcf {
    input:
      pav_conf = config,
      pav_asm = tar_asm.asm_tar,
      threads = "1",
      mem_gb = "16",
      sample = sample,
      contigInfo = data_ref_contig_table.contigInfo,
      finalBedOut = call_final_bed.bed
  }
  #call RunPaftools {
  #  input:
  #    hapOne = hapOne,
  #    hapTwo = hapTwo,
  #    ref = ref
  #}

  output {
    File pav_vcf = write_vcf.vcf
    #File paftools_vcf = RunPaftools.vcf
  }
}


task RunPaftools {
    input { 
        File hapOne
        File hapTwo
        File ref
    }
    Int disk_size_gb = ceil(size(hapOne, "GB") + size(hapTwo, "GB") + size(ref, "GB"))*10
    String docker_dir = "/assemblybased"
    String work_dir = "/cromwell_root/assemblybased"
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        MINIMAP2_ASM_FLAG="asm5"  # asm5/asm10/asm20 for ~0.1/1/5% sequence divergence

        cat ~{hapOne} ~{hapTwo} > h1_h2.fa
        ${TIME_COMMAND} minimap2 -t ${N_THREADS} -x ${MINIMAP2_ASM_FLAG} -c --cs ~{ref} h1_h2.fa -o h1_h2.paf
        rm -f h1_h2.fa
        ${TIME_COMMAND} sort --parallel=${N_THREADS} -k6,6 -k8,8n h1_h2.paf > h1_h2.sorted.paf
        rm -f h1_h2.paf
        ${TIME_COMMAND} paftools.js call -f ~{ref} h1_h2.sorted.paf > out.vcf
    >>>
    output {
        File vcf = work_dir + "/out.vcf"
    }
    runtime {
        docker: "fcunial/assemblybased"
        cpu: 32
        memory: "32GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 3
    }
}