### 1 : MODULE_1_DATA_ACQUISITION 
1. The ENA Throttling & Network Instability Problem

    The Error/Challenge: Initial data acquisition was extremely slow and prone to timeouts. The European Nucleotide Archive (ENA) aggressively throttles single-connection downloads, and standard network drops on a 100 Mbps Wi-Fi connection would cause the entire pipeline process to fail.

    Our Trial: We initially relied on standard, single-threaded download commands (like curl or wget) without robust error handling.

    The Solution: We implemented aria2c as our primary engine, configuring it to open 16 simultaneous connections per file (-x 16 -s 16). This successfully bypassed ENA's server-side throttling. We also parameterized this (aria2c_connections) in nextflow.config so it can be easily adjusted later.

2. The Missing Accession Problem (Single Point of Failure)

    The Error/Challenge: Relying strictly on the ENA FTP server meant that if an accession was delayed in mirroring from NCBI, or if the ENA server was temporarily down, the pipeline crashed.

    Our Trial: We attempted to handle missing files by just failing the pipeline and requiring manual intervention, which defeats the purpose of an automated pipeline.

    The Solution: We built a dual-strategy failover system. The script first attempts the ultra-fast ENA download. If (and only if) ENA fails or the file is missing, it automatically triggers a fallback to the NCBI SRA Toolkit (prefetch + fasterq-dump), safely zipping the output before moving on.

3. Hardware Saturation vs. Bottlenecking

    The Error/Challenge: Downloading massive .fastq.gz WGS files can easily overwhelm a local laptop’s SSD I/O and CPU, or conversely, fail to fully utilize the available bandwidth if run strictly sequentially.

    Our Trial: We initially hardcoded parallelization limits, which made the pipeline inflexible. We also attempted to count total reads (min_reads) during the download phase to filter out bad runs early.

    The Solution:

        Bandwidth Optimization: We introduced sra_max_forks (currently set to 4) to perfectly saturate the 100 Mbps network without melting the hardware.

        I/O Optimization: We completely removed the min_reads parameter from Module 1. Counting millions of reads required full decompression, causing a massive I/O bottleneck. We deferred this check to Module 2 (where fastp handles it natively), while keeping a lightweight Python script in Module 1 that only samples the first 2,000 reads to check average length (min_r2_len).

4. Silent Data Corruption

    The Error/Challenge: Network drops would sometimes result in partially downloaded FASTQ files. Nextflow would assume the process succeeded because a file was generated, passing corrupted data to downstream mapping tools and causing catastrophic failures later in the pipeline.

    Our Trial: Relying solely on the exit codes of the download tools.

    The Solution: We implemented hard integrity checks using gzip -t. Every single downloaded or merged file is tested for structural integrity. If a file is corrupt, it is immediately deleted, and Nextflow’s errorStrategy 'retry' automatically kicks in with an exponential back-off timer to try again.

Summary:

    "Module 1 is now a fault-tolerant data acquisition engine. Rather than relying on a single database or simple download commands, it actively routes around network failures, verifies data integrity at the byte level, and optimally balances our hardware resources. It guarantees that downstream tools will only ever receive uncorrupted, properly formatted FASTQ files."
### 2 : MODULE_3_MAP_READS
1) The Picard Memory Bottleneck & Duplicate Marking Optimization

    The Error/Challenge: Standard bioinformatics pipelines typically rely on Picard's MarkDuplicates for identifying PCR and optical duplicates. However, Picard is written in Java and is notoriously memory-intensive. On resource-constrained hardware (like a standard research laptop or Vivobook), it frequently crashes with Java heap-space out-of-memory errors (OutOfMemoryError), while also requiring heavy temporary file generation on disk.

    Our Trial: We initially looked into tuning Java virtual machine (JVM) memory allocation parameters (-Xmx) inside the configuration profile, but this either starved other processes of RAM or still failed when processing deeper genomic coverage files.

    The Solution: We completely eliminated Picard for duplicate marking and replaced it with a modern, pure-samtools compiled pipeline:

        samtools collate (groups reads by name efficiently)

        samtools fixmate (adds mate score tags required for duplicate flagging)

        samtools sort (re-sorts coordinates)

        samtools markdup -s (flags duplicates natively and outputs clean metrics)

    Why this is superior: Because samtools is written in low-level C, it runs significantly faster, has a tiny RAM footprint, bypasses Java memory crashes completely, and handles the operation natively without writing bloated intermediate files.