#!/usr/bin/env python
from __future__ import print_function

import re
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

from .utils import make_headers, Metric, exist_and_number

# Initialise the logger
import logging

log = logging.getLogger(__name__)


NAMESPACE = "CNV metrics"


class DragenCNVMetrics(BaseMultiqcModule):
    def add_cnv_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files("dragen/cnv"):
            s_name, data = parse_cnv_metrics_file(f)
            s_name = self.clean_s_name(s_name, f)
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, section="stats")
            data_by_sample[s_name] = data

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)
        if not data_by_sample:
            return set()

        # Write data to file
        self.write_data_file(data_by_sample, "dragen_cnv")

        all_metric_names = set()
        for sn, sdata in data_by_sample.items():
            for m in sdata.keys():
                all_metric_names.add(m)

        gen_stats_headers, cnv_table_headers = make_headers(all_metric_names, CNV_METRICS)

        self.general_stats_addcols(data_by_sample, gen_stats_headers, namespace=NAMESPACE)

        self.add_section(
            name="CNV metrics",
            anchor="dragen-cnv",
            description="""
            Summary of the CNV (Copy Number Variantion) data.
            """,
            plot=table.plot(data_by_sample, cnv_table_headers, pconfig={"namespace": NAMESPACE}),
        )
        return data_by_sample.keys()


CNV_METRICS = [
    Metric (
        m.id,
        m.title,
        in_genstats=m.in_genstats,
        in_own_tabl=m.in_own_tabl,
        descr=m.descr,
        unit=m.unit,
        namespace=m.namespace or NAMESPACE,
        the_higher_the_worse=m.the_higher_the_worse,
    )
    for m in [
        Metric(
            "Number of passing amplifications",
            "Number of passing amplifications",
            '#',
            '#',
            "",
            "Variants",
            "Number of passing amplifications",
            namespace="Variants",
        ),
        Metric(
            "Number of passing deletions",
            "Number of passing deletions",
            '#',
            '#',
            "",
            "Variants",
            "Number of passing deletions",
            namespace="Variants",
        ),
        Metric(
            "Bases in reference genome",
            "Bases in reference genome",
            None,
            'hid',
            "bases",
            "Bases in reference genome",
            namespace="bases"
        ),
        Metric(
            "Average alignment coverage over genome",
            "Average alignment coverage over genome",
            None,
            'hid',
            "",
            "Average alignment coverage over genome",
        ),
        Metric(
            "Number of alignment records",
            "Number of alignment records",
            None,
            'hid',
            "",
            "Number of alignment records",
        ),
        Metric(
            "Number of filtered records (total)",
            "Number of filtered records (total)",
            None,
            'hid',
            "",
            "Number of filtered records (total)",
        ),
        Metric(
            "Number of filtered records (duplicates)",
            "Number of filtered records (duplicates)",
            None,
            'hid',
            "",
            "Number of filtered records (duplicates)",
        ),
        Metric(
            "Number of filtered records (MAPQ)",
            "Number of filtered records (MAPQ)",
            None,
            'hid',
            "",
            "Number of filtered records (MAPQ)",
        ),
        Metric(
            "Number of filtered records (unmapped)",
            "Number of filtered records (unmapped)",
            None,
            'hid',
            "",
            "Number of filtered records (unmapped)",
        ),
        Metric(
            "Coverage MAD",
            "Coverage MAD",
            None,
            '#',
            "",
            "Coverage MAD",
        ),
        Metric(
            "Median Bin Count",
            "Median Bin Count",
            None,
            '#',
            "",
            "Median Bin Count",
        ),
        Metric(
            "Number of target intervals",
            "Number of target intervals",
            'hid',
            '#',
            "",
            "Number of target intervals",
        ),
        Metric(
            "Number of normal samples",
            "Number of normal samples",
            'hid',
            '#',
            "",
            "Number of normal samples",
        ),
        Metric(
            "Number of segments",
            "Number of segments",
            'hid',
            '#',
            "",
            "Number of segments",
        ),
        Metric(
            "Number of amplifications",
            "Number of amplifications",
            None,
            'hid',
            "Variants",
            "Number of amplifications",
            namespace="Variants"
        ),
        Metric(
            "Number of deletions",
            "Number of deletions",
            None,
            'hid',
            "",
            "Variants",
            "Number of deletions",
            namespace="Variants",
        ),
    ]

]


def parse_cnv_metrics_file(f):
    """
    Acrometrix_50ng_R2.cnv_metrics.csv
    SEX GENOTYPER,,Acrometrix_10ng_R2_JM,XY,0.8516
    SEX GENOTYPER,,<normal#>,XX,0.7738
    CNV SUMMARY,,Bases in reference genome,3613185058
    CNV SUMMARY,,Average alignment coverage over genome,12.94
    CNV SUMMARY,,Number of alignment records,326351771
    CNV SUMMARY,,Number of filtered records (total),7328181,2.25
    CNV SUMMARY,,Number of filtered records (duplicates),0,0.00
    CNV SUMMARY,,Number of filtered records (MAPQ),3665231,1.12
    CNV SUMMARY,,Number of filtered records (unmapped),3662950,1.12
    CNV SUMMARY,,Coverage MAD,0.12356
    CNV SUMMARY,,Median Bin Count,4.89
    CNV SUMMARY,,Number of target intervals,435022
    CNV SUMMARY,,Number of normal samples,29
    CNV SUMMARY,,Number of segments,12331
    CNV SUMMARY,,Number of amplifications,289
    CNV SUMMARY,,Number of deletions,171
    CNV SUMMARY,,Number of passing amplifications,35,12.11
    CNV SUMMARY,,Number of passing deletions,40,23.39
    """

    s_name = re.search(r"(.*)\.cnv_metrics.csv", f["fn"]).group(1)

    cnv_data = dict()
    genotyper_data = dict()

    for line in f["f"].splitlines():
        fields = line.split(",")
        analysis = fields[0]
        # sample = fields[1]
        metric = fields[2]
        value = fields[3]

        try:
            value = int(value)
        except ValueError:
            try:
                value = float(value)
            except ValueError:
                pass

        percentage = None
        if len(fields) > 4:  # percentage
            percentage = fields[4]
            try:
                percentage = float(percentage)
            except ValueError:
                pass

        if analysis == "SEX GENOTYPER":
            genotyper_data[metric] = value

        if analysis == "CNV SUMMARY":
            cnv_data[metric] = value
            if percentage is not None:
               cnv_data[metric + " pct"] = percentage


    data = cnv_data

    return s_name, data
