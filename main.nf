#!/usr/bin/env nextflow
/*
 * ont-amr-nf
 * Oxford Nanopore (R10.4.1) bacterial reads -> validated antimicrobial-resistance gene calls
 *
 * Design principle: the pipeline's job is not to produce a result, it is to produce a
 * result you can defend. Every isolate carries a role (test / control / titration) and
 * every AMR call passes a validation gate that can mark it uninterpretable.
 */

nextflow.enable.dsl = 2

include { FETCH_READS      } from './modules/fetch_reads'
include { NANOQ_FILTER     } from './modules/nanoq_filter'
include { MAKE_DECOY       } from './modules/make_decoy'
include { FLYE_ASSEMBLE    } from './modules/flye_assemble'
include { ASSEMBLY_STATS   } from './modules/assembly_stats'
include { AMRFINDERPLUS    } from './modules/amrfinderplus'
include { VALIDATION_GATE  } from './modules/validation_gate'
include { SUBSAMPLE_DEPTH  } from './modules/subsample_depth'
include { SHUFFLE_ASSEMBLY } from './modules/shuffle_assembly'
include { AGGREGATE_RESULTS } from './modules/aggregate'

// The caller-level control runs through the SAME process as every real isolate, only
// aliased so Nextflow can invoke it twice in one workflow. Same code, same parameters.
include { AMRFINDERPLUS as AMRFINDERPLUS_CONTROL } from './modules/amrfinderplus'

def helpMessage() {
    log.info """
    ont-amr-nf ${workflow.manifest.version}

    Usage:
      nextflow run . -profile conda --samplesheet assets/samplesheet.csv

    Key options:
      --samplesheet        CSV with sample_id,accession,species,role   [${params.samplesheet}]
      --outdir             output directory                            [${params.outdir}]
      --min_read_len       read length floor (bp)                      [${params.min_read_len}]
      --min_read_q         read quality floor                          [${params.min_read_q}]
      --min_depth_x        depth below which AMR absence is not called [${params.min_depth_x}]
      --make_decoy         generate the shuffled negative control      [${params.make_decoy}]
      --run_titration      run the depth-titration experiment          [${params.run_titration}]
      --titration_sample   sample_id to titrate                        [${params.titration_sample}]
      --titration_depths   comma-separated target depths               [${params.titration_depths}]

    Profiles: conda, docker, singularity, local_env, test
    """.stripIndent()
}

/*
 * Reads-to-calls for one channel of (meta, reads). Used for real isolates, for the
 * decoy control, and for every titration point, so all three are processed by
 * identical code — a control that took a different code path would prove nothing.
 */
workflow CALL_AMR {
    take:
    ch_reads

    main:
    NANOQ_FILTER(ch_reads)
    FLYE_ASSEMBLE(NANOQ_FILTER.out.reads)

    // Pair each assembly back with the filtered reads it came from, so depth can be
    // measured by remapping rather than assumed from input yield.
    ch_for_stats = FLYE_ASSEMBLE.out.assembly
        .join(NANOQ_FILTER.out.reads, by: 0)

    ASSEMBLY_STATS(ch_for_stats)
    AMRFINDERPLUS(FLYE_ASSEMBLE.out.assembly)

    // The gate needs the assembly status too: "no calls because nothing assembled" and
    // "no calls in a good assembly" are different findings and must not collapse.
    ch_gate = ASSEMBLY_STATS.out.stats
        .join(AMRFINDERPLUS.out.calls,   by: 0)
        .join(FLYE_ASSEMBLE.out.status,  by: 0)

    VALIDATION_GATE(ch_gate)

    emit:
    amr      = AMRFINDERPLUS.out.calls
    // The assembly as it was handed to the caller — the caller-level control shuffles
    // exactly this, so the control and the real call share an identical starting point.
    amr_input = FLYE_ASSEMBLE.out.assembly
    stats    = ASSEMBLY_STATS.out.stats
    verdict  = VALIDATION_GATE.out.verdict
    json     = VALIDATION_GATE.out.json
    versions = NANOQ_FILTER.out.versions
                 .mix(FLYE_ASSEMBLE.out.versions,
                      ASSEMBLY_STATS.out.versions,
                      AMRFINDERPLUS.out.versions)
}

workflow {

    if (params.help) {
        helpMessage()
        return
    }

    log.info """
    ================================================================
     ont-amr-nf ${workflow.manifest.version}
    ================================================================
     samplesheet   : ${params.samplesheet}
     outdir        : ${params.outdir}
     read filter   : >=${params.min_read_len} bp, >=Q${params.min_read_q}
     assembler     : flye ${params.flye_mode}
     depth floor   : ${params.min_depth_x}x (AMR absence below this is not reported)
     negative ctrl : ${params.make_decoy}
     titration     : ${params.run_titration}
    ================================================================
    """.stripIndent()

    // ---- samplesheet -> (meta, accession) ----
    ch_input = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .filter { row -> row.accession && row.accession != 'synthetic' }
        .map { row ->
            def meta = [
                id       : row.sample_id,
                accession: row.accession,
                species  : row.species,
                role     : row.role ?: 'test'
            ]
            tuple(meta, row.accession)
        }

    FETCH_READS(ch_input)
    ch_reads = FETCH_READS.out.reads

    // ---- negative control, derived from a real sample ----
    if (params.make_decoy) {
        ch_decoy_source = params.decoy_from
            ? ch_reads.filter { meta, r -> meta.id == params.decoy_from }
            : ch_reads.first()
        MAKE_DECOY(ch_decoy_source)
        // Relabel the decoy here rather than inside the process: its identity is a
        // workflow-level fact (role = negative_control), and keeping it out of the
        // script block lets the decoy share one code path with the real isolates.
        ch_decoy = MAKE_DECOY.out.reads
            .map { meta, reads ->
                def dmeta = [
                    id       : 'NEG_DECOY',
                    accession: 'synthetic',
                    species  : meta.species,
                    role     : 'negative_control',
                    shuffled_from: meta.id
                ]
                tuple(dmeta, reads)
            }
        ch_all = ch_reads.mix(ch_decoy)
    } else {
        ch_all = ch_reads
    }

    // ---- depth titration ----
    if (params.run_titration) {
        ch_tit_source = params.titration_sample
            ? ch_reads.filter { meta, r -> meta.id == params.titration_sample }
            : ch_reads.first()

        ch_depths = Channel.fromList(
            params.titration_depths.toString().split(',').collect { it.trim() as int } )
        ch_reps = Channel.fromList( (1..(params.titration_replicates as int)).toList() )

        ch_tit_jobs = ch_tit_source
            .combine(ch_depths)
            .combine(ch_reps)
            .map { meta, reads, depth, rep -> tuple(meta, reads, depth, rep) }

        SUBSAMPLE_DEPTH(ch_tit_jobs)
        // Relabel each titration point at the workflow level, same as the decoy, so the
        // process stays free of identity logic and every point joins the shared path.
        ch_tit = SUBSAMPLE_DEPTH.out.reads
            .map { meta, depth, rep, reads ->
                def tmeta = meta + [
                    id          : "${meta.id}_d${depth}_r${rep}",
                    parent_id   : meta.id,
                    target_depth: depth,
                    replicate   : rep,
                    role        : 'titration'
                ]
                tuple(tmeta, reads)
            }
        ch_all = ch_all.mix(ch_tit)
    }

    // ---- one code path for isolates, controls and titration points ----
    CALL_AMR(ch_all)

    // ---- caller-level control ----
    // The read-level decoy tests the ASSEMBLER: shuffled reads do not assemble, so the
    // caller only ever sees an empty FASTA and returns an empty call set trivially.
    // This control tests the CALLER: it shuffles bases within a real, good assembly —
    // preserving contig count, lengths and GC exactly — and sends the result through
    // the same AMRFinderPlus process with the same parameters. A call here is a call
    // driven by composition rather than gene identity.
    if (params.caller_control) {
        ch_cc_source = CALL_AMR.out.amr_input
            .filter { meta, asm -> meta.role == 'test' }
            .first()

        SHUFFLE_ASSEMBLY(ch_cc_source)

        ch_cc = SHUFFLE_ASSEMBLY.out.assembly
            .map { meta, asm ->
                def cmeta = [
                    id       : 'CALLER_CONTROL',
                    accession: 'synthetic',
                    species  : meta.species,
                    role     : 'caller_control',
                    shuffled_from: meta.id
                ]
                tuple(cmeta, asm)
            }

        AMRFINDERPLUS_CONTROL(ch_cc)
        ch_cc_calls = AMRFINDERPLUS_CONTROL.out.calls
    } else {
        ch_cc_calls = Channel.empty()
    }

    AGGREGATE_RESULTS(
        CALL_AMR.out.amr.map    { meta, f -> f }
            .mix(ch_cc_calls.map { meta, f -> f })
            .collect(),
        CALL_AMR.out.stats.map  { meta, f -> f }.collect(),
        CALL_AMR.out.json.collect()
    )

    CALL_AMR.out.versions
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info",
                     sort: true, keepHeader: false)

    workflow.onComplete = { completionMessage() }
}

// Nextflow 26.x requires the completion handler to be declared inside a workflow
// block rather than as a top-level statement.
def completionMessage() {
    log.info """
    ================================================================
     completed : ${workflow.success ? 'OK' : 'FAILED'}
     duration  : ${workflow.duration}
     results   : ${params.outdir}/
       - amr_calls.tsv           every resistance determinant, per sample
       - validation_summary.tsv  PASS/FAIL per sample with the failing checks
       - run_summary.md          human-readable overview
       - pipeline_info/          execution report, trace, DAG, tool versions
    ================================================================
    """.stripIndent()
}
