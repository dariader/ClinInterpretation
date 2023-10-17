#!/usr/bin/env python

# Кушает хромосому в формате ChrNN, координату в 38-й сборке, референсный и альтернативный аллель
def gnomad_checker(chrom, coord, ref, alt):


    chrom_n = chrom[3:]

    url = "https://gnomad.broadinstitute.org/api/"

    headers = {"Content-Type": "application/json",
               "Host": "gnomad.broadinstitute.org",
               "Origin": "https://gnomad.broadinstitute.org",
               "Referer": f"https://gnomad.broadinstitute.org/variant/{chrom_n}-{coord}-{ref}-{alt}?dataset=gnomad_r3"}


    json_data = {"operationName":"GnomadVariant",
                 "query":"\nquery GnomadVariant($variantId: String!, $datasetId: DatasetId!, $referenceGenome: ReferenceGenomeId!, $includeLocalAncestry: Boolean!, $includeLiftoverAsSource: Boolean!, $includeLiftoverAsTarget: Boolean!) {\n  variant(variantId: $variantId, dataset: $datasetId) {\n    variant_id\n    reference_genome\n    chrom\n    pos\n    ref\n    alt\n    caid\n    colocated_variants\n    coverage {\n      exome {\n        mean\n      }\n      genome {\n        mean\n      }\n    }\n    multi_nucleotide_variants {\n      combined_variant_id\n      changes_amino_acids\n      n_individuals\n      other_constituent_snvs\n    }\n    exome {\n      ac\n      an\n      ac_hemi\n      ac_hom\n      faf95 {\n        popmax\n        popmax_population\n      }\n      filters\n      populations {\n        id\n        ac\n        an\n        ac_hemi\n        ac_hom\n      }\n      local_ancestry_populations @include(if: $includeLocalAncestry) {\n        id\n        ac\n        an\n      }\n      age_distribution {\n        het {\n          bin_edges\n          bin_freq\n          n_smaller\n          n_larger\n        }\n        hom {\n          bin_edges\n          bin_freq\n          n_smaller\n          n_larger\n        }\n      }\n      quality_metrics {\n        allele_balance {\n          alt {\n            bin_edges\n            bin_freq\n            n_smaller\n            n_larger\n          }\n        }\n        genotype_depth {\n          all {\n            bin_edges\n            bin_freq\n            n_smaller\n            n_larger\n          }\n          alt {\n            bin_edges\n            bin_freq\n            n_smaller\n            n_larger\n          }\n        }\n        genotype_quality {\n          all {\n            bin_edges\n            bin_freq\n            n_smaller\n            n_larger\n          }\n          alt {\n            bin_edges\n            bin_freq\n            n_smaller\n            n_larger\n          }\n        }\n        site_quality_metrics {\n          metric\n          value\n        }\n      }\n    }\n    genome {\n      ac\n      an\n      ac_hemi\n      ac_hom\n      faf95 {\n        popmax\n        popmax_population\n      }\n      filters\n      populations {\n        id\n        ac\n        an\n        ac_hemi\n        ac_hom\n      }\n      local_ancestry_populations @include(if: $includeLocalAncestry) {\n        id\n        ac\n        an\n      }\n      age_distribution {\n        het {\n          bin_edges\n          bin_freq\n          n_smaller\n          n_larger\n        }\n        hom {\n          bin_edges\n          bin_freq\n          n_smaller\n          n_larger\n        }\n      }\n      quality_metrics {\n        allele_balance {\n          alt {\n            bin_edges\n            bin_freq\n            n_smaller\n            n_larger\n          }\n        }\n        genotype_depth {\n          all {\n            bin_edges\n            bin_freq\n            n_smaller\n            n_larger\n          }\n          alt {\n            bin_edges\n            bin_freq\n            n_smaller\n            n_larger\n          }\n        }\n        genotype_quality {\n          all {\n            bin_edges\n            bin_freq\n            n_smaller\n            n_larger\n          }\n          alt {\n            bin_edges\n            bin_freq\n            n_smaller\n            n_larger\n          }\n        }\n        site_quality_metrics {\n          metric\n          value\n        }\n      }\n    }\n    flags\n    lof_curations {\n      gene_id\n      gene_symbol\n      verdict\n      flags\n      project\n    }\n    rsids\n    transcript_consequences {\n      domains\n      gene_id\n      gene_version\n      gene_symbol\n      hgvs\n      hgvsc\n      hgvsp\n      is_canonical\n      is_mane_select\n      is_mane_select_version\n      lof\n      lof_flags\n      lof_filter\n      major_consequence\n      polyphen_prediction\n      sift_prediction\n      transcript_id\n      transcript_version\n    }\n    in_silico_predictors {\n      id\n      value\n      flags\n    }\n    non_coding_constraint {\n      start\n      stop\n      possible\n      observed\n      expected\n      oe\n      z\n    }\n  }\n\n  clinvar_variant(variant_id: $variantId, reference_genome: $referenceGenome) {\n    clinical_significance\n    clinvar_variation_id\n    gold_stars\n    last_evaluated\n    review_status\n    submissions {\n      clinical_significance\n      conditions {\n        name\n        medgen_id\n      }\n      last_evaluated\n      review_status\n      submitter_name\n    }\n  }\n\n  liftover(source_variant_id: $variantId, reference_genome: $referenceGenome) @include(if: $includeLiftoverAsSource) {\n    liftover {\n      variant_id\n      reference_genome\n    }\n    datasets\n  }\n\n  liftover_sources: liftover(liftover_variant_id: $variantId, reference_genome: $referenceGenome) @include(if: $includeLiftoverAsTarget) {\n    source {\n      variant_id\n      reference_genome\n    }\n    datasets\n  }\n\n  meta {\n    clinvar_release_date\n  }\n}\n",
                 "variables":{"datasetId":"gnomad_r3",
                              "includeLocalAncestry":True,
                              "includeLiftoverAsSource":False,
                              "includeLiftoverAsTarget":True,
                              "referenceGenome":"GRCh38",
                              "variantId": f"{chrom_n}-{coord}-{ref}-{alt}"}}


    gnomad_response = requests.post(url=url, headers=headers, json=json_data)
    res = gnomad_response.json()
    if res["data"]["variant"] is None:
        af = 0
    else:
        gnomad_female = res["data"]["variant"]["genome"]["populations"][30]
        gnomad_male = res["data"]["variant"]["genome"]["populations"][31]

        gnomADv3_AC_XX = gnomad_female["ac"]
        gnomADv3_AC_XY = gnomad_male["ac"]


    return gnomADv3_AC_XX, gnomADv3_AC_XY
