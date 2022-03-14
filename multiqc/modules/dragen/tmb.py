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


NAMESPACE = "TMB metrics"


class DragenTMBMetrics(BaseMultiqcModule):
    def add_tmb_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files("dragen/tmb"):
            s_name, data = parse_tmb_metrics_file(f)
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
        self.write_data_file(data_by_sample, "dragen_tmb")

        all_metric_names = set()
        for sn, sdata in data_by_sample.items():
            for m in sdata.keys():
                all_metric_names.add(m)

        gen_stats_headers, tmb_table_headers = make_headers(all_metric_names, TMB_METRICS)

        self.general_stats_addcols(data_by_sample, gen_stats_headers, namespace=NAMESPACE)

        self.add_section(
            name="TMB",
            anchor="dragen-tmb",
            description="""
            Estimation of the Tumor Mutational Burden.
            """,
            plot=table.plot(data_by_sample, tmb_table_headers, pconfig={"namespace": NAMESPACE}),
        )
        return data_by_sample.keys()


TMB_METRICS = [
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
            "Total Input Variant Count",
            "Total Input Variant Count",
            None,
            'hid',
            "Variants",
            "Total Input Variant Count",
            namespace="Variants"
        ),
        Metric(
            "Total Input Variant Count in TMB region",
            "Total Input Variant Count in TMB region",
            'None',
            '#',
            "Variants",
            "Total Input Variant Count in TMB region",
            namespace="Variants"
        ),
        Metric(
            "Filtered Variant Count",
            "Filtered Variant Count",
            None,
            '#',
            "Variants",
            "Filtered Variant Count",
            namespace="Variants"
        ),
        Metric(
            "Filtered Nonsyn Variant Count",
            "Filtered Nonsyn Variant Count",
            None,
            '#',
            "Variants",
            "Filtered Nonsyn Variant Count",
            namespace="Variants"
        ),
        Metric(
            "Eligible Region (MB)",
            "Eligible Region (MB)",
            None,
            'hid',
            "",
            "Eligible Region (MB)",
        ),
        Metric(
            "TMB",
            "TMB",
            'hid',
            '#',
            "",
            "Tumor Mutational Burden"
        ),
        Metric(
            "Nonsyn TMB",
            "Nonsyn TMB",
            '#',
            '#',
            "",
            "Nonsynonymous Tumor Mutational Burden"
        ),
    ]

]


def parse_tmb_metrics_file(f):
    """
    Acrometrix_50ng_R1.tmb.metrics.csv

    TMB SUMMARY,,Total Input Variant Count,26546
    TMB SUMMARY,,Total Input Variant Count in TMB region,26311
    TMB SUMMARY,,Filtered Variant Count,456
    TMB SUMMARY,,Filtered Nonsyn Variant Count,393
    TMB SUMMARY,,Eligible Region (MB),35.46
    TMB SUMMARY,,TMB,12.86
    TMB SUMMARY,,Nonsyn TMB,11.08
    """

    s_name = re.search(r"(.*)\.tmb.metrics.csv", f["fn"]).group(1)

    data = defaultdict(dict)

    for line in f["f"].splitlines():
        _, _, metric, stat = line.split(",")
        try:
            stat = float(stat)
        except ValueError:
            pass
        data[metric] = stat

    return s_name, data
