process FETCH_READS {
    tag   "${meta.id}"
    label 'process_low'
    conda "conda-forge::python=3.11"
    publishDir "${params.outdir}/reads", mode: params.publish_mode, pattern: "*.md5"
    maxForks 2          // be polite to NCBI
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: reads
    path "${meta.id}.md5",                        emit: checksum

    script:
    """
    #!/usr/bin/env python3
    import urllib.request, hashlib, os, shutil, sys, time

    # If a local copy already exists, use it: re-runs stay offline, the checksum record
    # stays comparable across runs, and NCBI is not re-hit for data already on disk.
    # --reads_dir overrides the default location (e.g. an external volume).
    cands = [
        os.path.join("${params.reads_dir ?: ''}", "${meta.id}.fastq.gz"),
        os.path.join("${projectDir}", "data", "raw", "${meta.id}.fastq.gz"),
    ]
    local = next((p for p in cands if p and os.path.exists(p)), None)
    dest  = "${meta.id}.fastq.gz"

    if local and os.path.getsize(local) > 1_000_000:
        # Copy rather than symlink: the reads may live on a filesystem that does not
        # support symlinks (e.g. a FAT-formatted external volume).
        shutil.copyfile(local, dest)
        src = "local:" + local
    else:
        url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=${accession}"
        src = url
        last = None
        for attempt in range(4):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "ont-amr-nf"})
                with urllib.request.urlopen(req, timeout=1800) as r, open(dest, "wb") as out:
                    while True:
                        chunk = r.read(1 << 20)
                        if not chunk:
                            break
                        out.write(chunk)
                break
            except Exception as e:
                last = e
                time.sleep(10 * (attempt + 1))
        else:
            sys.exit(f"failed to fetch ${accession}: {last}")

    if os.path.getsize(dest) < 1_000_000:
        sys.exit("fetched file is implausibly small for an ONT WGS run")

    h = hashlib.md5()
    with open(dest, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)

    with open("${meta.id}.md5", "w") as fh:
        fh.write(f"{h.hexdigest()}\\t${meta.id}\\t${accession}\\t{src}\\n")
    """

    stub:
    """
    echo -e "@r1\\nACGT\\n+\\nIIII" | gzip > ${meta.id}.fastq.gz
    touch ${meta.id}.md5
    """
}
