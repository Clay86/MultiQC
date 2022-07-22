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


NAMESPACE = "time metrics"


class DragenTimeMetrics(BaseMultiqcModule):
    def add_time_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files("dragen/time"):
            s_name, data = parse_time_metrics_file(f)
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
        self.write_data_file(data_by_sample, "dragen_time")

        all_metric_names = set()
        for sn, sdata in data_by_sample.items():
            for m in sdata.keys():
                all_metric_names.add(m)

        gen_stats_headers, time_table_headers = make_headers(all_metric_names, TIME_METRICS)

        for s in data_by_sample:
            for m in all_metric_names:
                if m not in data_by_sample[s].keys():
                    data_by_sample[s][m] = ""

        self.general_stats_addcols(data_by_sample, gen_stats_headers, namespace=NAMESPACE)


        self.add_section(
            name="Time Metrics",
            anchor="dragen-time",
            description="""
            Duration of the pipeline exectution.
            """,
            plot=table.plot(data_by_sample, time_table_headers, pconfig={"namespace": NAMESPACE}),
        )
        return data_by_sample.keys()


TIME_METRICS = [
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
            "Total runtime",
            "Total runtime",
            'hid',
            '#',
            descr="Total runtime"
        ),
        Metric(
            "Time loading reference",
            "Time loading reference",
            None,
            'hid',
            descr="Time loading reference"
        ),
        Metric(
            "Time aligning reads",
            "Time aligning reads",
            None,
            'hid',
            descr="Time aligning reads"
        ),
        Metric(
            "Time duplicate marking",
            "Time duplicate marking",
            None,
            'hid',
            descr="Time duplicate marking"
        ),
        Metric(
            "Time sorting and marking duplicates",
            "Time sorting and marking duplicates",
            None,
            'hid',
            descr="Time sorting and marking duplicates"
        ),
        Metric(
            "Time DRAGStr calibration",
            "Time DRAGStr calibration",
            None,
            'hid',
            descr="Time DRAGStr calibration"
        ),
        Metric(
            "Time saving map/align output",
            "Time saving map/align output",
            None,
            'hid',
            descr="Time saving map/align output"
        ),
        Metric(
            "Time partial reconfiguration",
            "Time partial reconfiguration",
            None,
            'hid',
            descr="Time partial reconfiguration"
        ),
        Metric(
            "Time processing depth of coverage",
            "Time processing depth of coverage",
            None,
            'hid',
            descr="Time processing depth of coverage"
        ),
        Metric(
            "Time variant calling",
            "Time variant calling",
            'hid',
            '#',
            descr="Time variant calling"
        ),
        Metric(
            "Time calculating target counts",
            "Time calculating target counts",
            'hid',
            '#',
            descr="Time calculating target counts"
        ),
        Metric(
            "Time correcting GC bias",
            "Time correcting GC bias",
            'hid',
            '#',
            descr="Time correcting GC bias"
        ),
        Metric(
            "Time normalizing case sample",
            "Time normalizing case sample",
            None,
            'hid',
            descr="Time normalizing case sample"
        ),
        Metric(
            "Time performing segmentation",
            "Time performing segmentation",
            'hid',
            '#',
            descr="Time performing segmentation"
        ),
        Metric(
            "Time generating CNV calls",
            "Time generating CNV calls",
            'hid',
            '#',
            descr="Time generating CNV calls"
        ),
        Metric(
            "Time generating CNV track files",
            "Time generating CNV track files",
            None,
            'hid',
            descr="Time generating CNV track files"
        ),
        Metric(
            "Time partitioning",
            "Time partitioning",
            None,
            'hid',
            descr="Time partitioning"
        ),
        Metric(
            "Time annotating outputs",
            "Time annotating outputs",
            'hid',
            '#',
            descr="Time annotating outputs"
        ),
        Metric(
            "Time calculating TMB",
            "Time calculating TMB",
            'hid',
            '#',
            descr="Time calculating TMB"
        ),
    ]

]


def parse_time_metrics_file(f):
    """
    """

    s_name = re.search(r"(.*)\.time.metrics.csv", f["fn"]).group(1)

    data = defaultdict(dict)

    for line in f["f"].splitlines():
        _, _, metric, _, stat = line.split(",")
        try:
            stat = float(stat)
        except ValueError:
            pass
        data[metric] = stat

    return s_name, data
