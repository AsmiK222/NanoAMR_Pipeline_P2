#!/usr/bin/env bash
set -euo pipefail

ISOLATES=()
for i in $(seq -w 1 23); do
    ISOLATES+=("SA_${i}")
done

REF_ACCESSION="GCF_000013425.1"
REF_ASSEMBLY_NAME="ASM1342v1"
REF_DIR="reference"
ALIGN_DIR="alignments"
QC_DIR="alignment_qc"
CONSENSUS_DIR="consensus"
ASSEMBLY_DIR="assembly"

GENERATE_CONSENSUS=true
RUN_FLYE=false

MAMBA_BIN="$HOME/tools/bin/micromamba"

mkdir -p "$REF_DIR" "$ALIGN_DIR" "$QC_DIR"
[ "$GENERATE_CONSENSUS" = true ] && mkdir -p "$CONSENSUS_DIR"
[ "$RUN_FLYE" = true ] && mkdir -p "$ASSEMBLY_DIR"

REF_FASTA="$REF_DIR/S_aureus_ref.fasta"
if [ ! -f "$REF_FASTA" ]; then
    echo ">>> Downloading reference genome ($REF_ACCESSION)..."
    FTP_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/${REF_ACCESSION:0:3}/${REF_ACCESSION:4:3}/${REF_ACCESSION:7:3}/${REF_ACCESSION:10:3}/${REF_ACCESSION}_${REF_ASSEMBLY_NAME}/${REF_ACCESSION}_${REF_ASSEMBLY_NAME}_genomic.fna.gz"
    curl -f -o "$REF_FASTA.gz" "$FTP_URL"
    gunzip "$REF_FASTA.gz"
    samtools faidx "$REF_FASTA"
else
    echo ">>> Reference already present at $REF_FASTA, skipping download."
    [ -f "$REF_FASTA.fai" ] || samtools faidx "$REF_FASTA"
fi

# --- Consensus command: use the modern micromamba samtools (fast, built for
# long-read data) if available. This is what fixed the 43-minute bcftools
# hang -- do not fall back to bcftools unless this genuinely isn't installed.
samtools_new() {
    "$MAMBA_BIN" run -n bioenv samtools "$@"
}

CONSENSUS_METHOD=""
if [ "$GENERATE_CONSENSUS" = true ]; then
    if [ -x "$MAMBA_BIN" ]; then
        CONSENSUS_METHOD="samtools_new"
        echo ">>> Consensus method: modern samtools consensus (via micromamba) — fast path"
    elif samtools 2>&1 | grep -qi "consensus"; then
        CONSENSUS_METHOD="samtools_system"
        echo ">>> Consensus method: system samtools consensus"
    else
        echo "!!! No fast consensus method available. Falling back to bcftools —"
        echo "!!! WARNING: this took 43+ minutes per isolate last time on nanopore data. Not recommended."
        CONSENSUS_METHOD="bcftools"
    fi
fi

TOTAL=${#ISOLATES[@]}
COUNT=0
for ISOLATE in "${ISOLATES[@]}"; do
    COUNT=$((COUNT + 1))
    FASTQ="raw_data/${ISOLATE}/${ISOLATE}_reads.fastq.gz"
    BAM="$ALIGN_DIR/${ISOLATE}.sorted.bam"
    QC_FILE="$QC_DIR/${ISOLATE}_alignment_qc.txt"
    CONSENSUS_FASTA="$CONSENSUS_DIR/${ISOLATE}_consensus.fasta"

    echo ""
    echo "=== [$COUNT/$TOTAL] $ISOLATE ==="

    if [ ! -f "$FASTQ" ]; then
        echo "!!! Missing $FASTQ — skipping $ISOLATE"
        continue
    fi

    if [ -f "$BAM" ] && [ -f "$QC_FILE" ]; then
        echo ">>> Alignment already done, skipping minimap2/sort/QC."
    else
        SECONDS=0
        SAM="$ALIGN_DIR/${ISOLATE}.sam"
        echo ">>> Aligning with minimap2 (map-ont preset)..."
        minimap2 -ax map-ont "$REF_FASTA" "$FASTQ" > "$SAM" 2> "$ALIGN_DIR/${ISOLATE}.minimap2.log"
        echo ">>> Sorting and indexing..."
        samtools sort -@ 4 -o "$BAM" "$SAM"
        samtools index "$BAM"
        rm -f "$SAM"
        echo ">>> Alignment QC..."
        {
            echo "=== $ISOLATE alignment QC ==="
            samtools flagstat "$BAM"
            echo "--- mean depth ---"
            samtools depth -a "$BAM" | awk '{sum+=$3; n++} END {if (n>0) print sum/n, "x mean coverage"}'
        } > "$QC_FILE"
        cat "$QC_FILE"
        echo ">>> Alignment done in ${SECONDS}s"
    fi

    if [ "$GENERATE_CONSENSUS" = true ]; then
        if [ -f "$CONSENSUS_FASTA" ]; then
            echo ">>> Consensus already exists, skipping."
        else
            SECONDS=0
            echo ">>> Generating consensus FASTA (method: $CONSENSUS_METHOD)..."
            case "$CONSENSUS_METHOD" in
                samtools_new)
                    samtools_new consensus "$BAM" -o "$CONSENSUS_FASTA"
                    ;;
                samtools_system)
                    samtools consensus "$BAM" -o "$CONSENSUS_FASTA"
                    ;;
                bcftools)
                    VCF="$CONSENSUS_DIR/${ISOLATE}_calls.vcf.gz"
                    bcftools mpileup -Ou --max-depth 200 --ff SECONDARY,SUPPLEMENTARY -q 20 -f "$REF_FASTA" "$BAM" 2>/dev/null | bcftools call -mv -Oz --threads 4 -o "$VCF"
                    bcftools index -f "$VCF"
                    bcftools consensus -f "$REF_FASTA" "$VCF" > "$CONSENSUS_FASTA"
                    ;;
            esac
            echo ">>> Consensus done in ${SECONDS}s"
        fi
    fi

    if [ "$RUN_FLYE" = true ]; then
        echo ">>> Running Flye de novo assembly..."
        flye --nano-raw "$FASTQ" --out-dir "$ASSEMBLY_DIR/${ISOLATE}" --threads 4
    fi
done

echo ""
echo ">>> Stage 2 complete."
