process CHANGEO_PARSEDB {
    tag "$meta.id"
    label 'process_low'
    label 'changeo'

    input:
    tuple val(meta), path(db)   // AIRR TSV from MakeDb

    output:
    tuple val(meta), path("*_parse-pass.tsv"), emit: tab
    path "versions.yml",                       emit: versions

    script:
    """
    # 1. Filter productive sequences only (Briney uses productive-only)
    ParseDb.py select \\
        -d ${db} \\
        -f productive \\
        -u T \\
        --outname ${meta.id}-productive

    # 2. Strip allele from v_call / j_call into new gene-level columns
    #    Briney clonotype = (V_gene, J_gene, CDR3_aa) — gene-level, not allele
    python3 - "${meta.id}-productive_parse-select.tsv" "${meta.id}_parse-pass.tsv" <<'PYEOF'
import sys
import pandas as pd

in_tsv, out_tsv = sys.argv[1], sys.argv[2]
df = pd.read_csv(in_tsv, sep='\t')

def gene_only(call):
    if pd.isna(call):
        return call
    # IgBLAST may emit comma-separated ambiguous calls (e.g. "IGHV1-2*02,IGHV1-2*04")
    # Strip "*allele" from each, then collapse duplicates
    genes = [g.split('*')[0] for g in str(call).split(',')]
    return ','.join(sorted(set(genes)))

df['v_call_gene'] = df['v_call'].map(gene_only)
df['j_call_gene'] = df['j_call'].map(gene_only)
df.to_csv(out_tsv, sep='\t', index=False)
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        changeo: \$(python3 -c "import changeo; print(changeo.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_parse-pass.tsv
    echo '"${task.process}":' > versions.yml
    echo '    changeo: 1.3.0' >> versions.yml
    """
}
