import logging

from .coverage_hist import DragenCoverageHist
from .coverage_metrics import DragenCoverageMetrics
from .coverage_per_contig import DragenCoveragePerContig
from .fragment_length import DragenFragmentLength
from .mapping_metrics import DragenMappingMetrics
from .ploidy_estimation_metrics import DragenPloidyEstimationMetrics
from .vc_metrics import DragenVCMetrics
from .tmb import DragenTMBMetrics
from .coverage_metrics_covRegion import DragenCoverageRegionMetrics
from .cnv import DragenCNVMetrics
from .msi import DragenMsiMetrics
from .gc_metrics import DragenGcMetrics
from .rna_quant_metrics import DragenRnaQuantMetrics
from .rna_transcript_cov import DragenRnaTranscriptCoverage
from .sc_atac_metrics import DragenScAtacMetrics
from .sc_rna_metrics import DragenScRnaMetrics
from .time_metrics import DragenTimeMetrics
from .trimmer_metrics import DragenTrimmerMetrics

log = logging.getLogger(__name__)


class MultiqcModule(
    DragenMappingMetrics,
    DragenFragmentLength,
    DragenPloidyEstimationMetrics,
    DragenVCMetrics,
    DragenCoveragePerContig,
    DragenCoverageMetrics,
    DragenCoverageHist,
    DragenTMBMetrics,
    DragenCoverageRegionMetrics,
    DragenCNVMetrics,
    DragenMsiMetrics,
    DragenGcMetrics,
    DragenRnaQuantMetrics,
    DragenRnaTranscriptCoverage,
    DragenScAtacMetrics,
    DragenScRnaMetrics,
    DragenTimeMetrics,
    DragenTrimmerMetrics

):
    """DRAGEN provides a number of different pipelines and outputs, including base calling, DNA and RNA alignment,
    post-alignment processing and variant calling, covering virtually all stages of typical NGS data processing.
    However, it can be treated as a fast aligner with additional features on top, as users will unlikely use any
    features without enabling DRAGEN mapping. So we will treat this module as an alignment tool module and
    place it accordingly in the module_order list, in docs, etc.

    The QC metrics DRAGEN generates resemble those of samtools-stats, qualimap, mosdepth, bcftools-stats and alike.
    Whenever possible, the visual output is made similar to those modules.

    Note that this MultiQC module supports some of DRAGEN output but not all. Contributions are welcome!

    The code is structured in a way so every mix-in parses one type of QC file that DRAGEN generates
    (e.g. *.mapping_metrics.csv, *.wgs_fine_hist_normal.csv, etc). If a corresponding file is found, a mix-in adds
    a section into the report.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="DRAGEN",
            anchor="DRAGEN",
            target="DRAGEN",
            href="https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html",
            info=""" is a Bio-IT Platform that provides ultra-rapid secondary analysis of sequencing data
                     using field-programmable gate array technology (FPGA).""",
            # Can't find a DOI // doi=
        )

        samples_found = set()

        # <output prefix>.(wgs|target_bed)_coverage_metrics_(tumor|normal)?.csv
        # general stats table and a dedicated table
        samples_found |= self.add_coverage_region_metrics()

        # <output prefix>.(wgs|target_bed)_fine_hist_(tumor|normal)?.csv
        # coverage distribution and cumulative coverage plots
        samples_found |= self.add_coverage_hist()

        # <output prefix>.(wgs|target_bed)_coverage_metrics_(tumor|normal)?.csv
        # general stats table and a dedicated table
        samples_found |= self.add_coverage_metrics()

        # <output prefix>.(wgs|target_bed)_contig_mean_cov_(tumor|normal)?.csv
        # a histogram like in mosdepth, with each chrom as a category on X axis, plus a category
        # for autosomal chromosomes average
        samples_found |= self.add_coverage_per_contig()

        # general stats table, a dedicated table, and a few barplots
        # <output prefix>.mapping_metrics.csv
        samples_found |= self.add_mapping_metrics()

        # a histogram plot
        # <output prefix>.fragment_length_hist.csv
        samples_found |= self.add_fragment_length_hist()
        
        # <output prefix>.ploidy_estimation_metrics.csv
        # a "Ploidy estimation" column in the general stats table
        samples_found |= self.add_ploidy_estimation_metrics()

        # <output prefix>.vc_metrics.csv
        # a dedicated table and the total number of Variants into the general stats table
        samples_found |= self.add_vc_metrics()

        # CNV estimation
        # <output prefix>.cnv.metrics.csv
        samples_found |= self.add_cnv_metrics()

        # TMB estimation
        # <output prefix>.tmb.metrics.csv
        samples_found |= self.add_tmb_metrics()

        # MSI estimation
        # <output prefix>.microsat_output.json
        samples_found |= self.add_msi_metrics()

        # gc metrics
        # <output prefix>.gc_metrics.csv
        samples_found |= self.add_gc_metrics_hist()

        # rna quant metrics
        # <output prefix>.quant.metrics.csv
        samples_found |= self.add_rna_metrics()

        # rna transcript coverage
        # <output prefix>.quant.transcript_coverage.txt
        samples_found |= self.add_rna_transcript_coverage()

        # single cell atac metrics
        # <output prefix>.scATAC.metrics.csv
        samples_found |= self.add_sc_atac_metrics()

        # single cell rna metrics
        # <output prefix>.scRNA*metrics.csv
        samples_found |= self.add_sc_rna_metrics()

        # dragen time metrics
        # <output prefix>.time_metrics.csv
        samples_found |= self.add_time_metrics()

        #dragen trimmer metrics
        # <output prefix>.trimmer_metrics.csv
        samples_found |= self.add_trimmer_metrics()

        if len(samples_found) == 0:
            raise UserWarning
        log.info("Found {} reports".format(len(samples_found)))
