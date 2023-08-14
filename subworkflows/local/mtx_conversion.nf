/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { MTX_TO_H5AD   }             from '../../modules/local/mtx_to_h5ad.nf'
include { CONCAT_H5AD   }             from '../../modules/local/concat_h5ad.nf'
include { MTX_TO_SEURAT }             from '../../modules/local/mtx_to_seurat.nf'

workflow MTX_CONVERSION {

    take:
    mtx_matrices
    samplesheet
    txp2gene
    star_index

    main:
        ch_versions = Channel.empty()

        // Prepare Seurat mtx input, we let pass the files in filtered.../ *and* the filtered....h5
        mych_matrices = mtx_matrices.map { meta, ofiles ->
                                           def ffiles = ofiles.findAll { it =~ /filtered_feature_bc_matrix/ ||
                                                                         it =~ /count\/counts_unfiltered\// ||
                                                                         it =~ /af_quant\/alevin\//         ||
                                                                         it =~ /\.Solo\.out\/Gene/ }
                                           [meta, ffiles] }
        mych_matrices.dump(tag: 'MTX_CONVERSION(mod):PreH5AD', pretty: true)

        //
        // Convert matrix to h5ad
        //
        MTX_TO_H5AD (
            mych_matrices,
            txp2gene,
            star_index
        )

        //
        // Concat sample-specific h5ad in one
        //
        CONCAT_H5AD (
            MTX_TO_H5AD.out.h5ad.collect(), // gather all sample-specific files
            samplesheet
        )

        //
        // Convert matrix do seurat
        //
        MTX_TO_SEURAT (
            mych_matrices
        )

        //TODO CONCAT h5ad and MTX to h5ad should also have versions.yaml output
        ch_version = ch_versions.mix(MTX_TO_H5AD.out.versions, MTX_TO_SEURAT.out.versions)

    emit:
    ch_versions
    counts = MTX_TO_H5AD.out.counts

}
