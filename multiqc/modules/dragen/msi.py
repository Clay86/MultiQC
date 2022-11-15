#!/usr/bin/env python
from __future__ import print_function

import re
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table
import json

from .utils import make_headers, Metric, exist_and_number

# Initialise the logger
import logging

log = logging.getLogger(__name__)


NAMESPACE = "msi metrics"


class DragenMsiMetrics(BaseMultiqcModule):
    def add_msi_metrics(self):
        data_by_sample = dict()
        
        for f in self.find_log_files("dragen/msi"):
            s_name, data = parse_msi_metrics_file(f)
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
        self.write_data_file(data_by_sample, "dragen_msi")

        all_metric_names = set()
        for sn, sdata in data_by_sample.items():
            for m in sdata.keys():
                all_metric_names.add(m)

        gen_stats_headers, msi_table_headers = make_headers(all_metric_names, msi_METRICS)

        self.general_stats_addcols(data_by_sample, gen_stats_headers, namespace=NAMESPACE)

        self.add_section(
            name="MSI",
            anchor="dragen-msi",
            description="""
            Estimation of the Microsatellite Stability.
            """,
            plot=table.plot(data_by_sample, msi_table_headers, pconfig={"namespace": NAMESPACE}),
        )
        return data_by_sample.keys()


msi_METRICS = [
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
            "TotalMicrosatelliteSitesAssessed",
            "Total Microsatellite Sites Assessed",
            'hid',
            '#',
            descr="Total Microsatellite Sites Assessed"
        ),
        Metric(
            "TotalMicrosatelliteSitesUnstable",
            "Total Microsatellite Sites Unstable",
            '#',
            '#',
            descr="Total Microsatellite Sites Unstable"
        ),
        Metric(
            "PercentageUnstableSites",
            "Percentage Unstable MSI Sites",
            'hid',
            '#',
            "%",
            descr="Percentage Unstable MSI Sites",
        ),
        Metric(
            "ResultIsValid",
            "Result Is Valid",
            None,
            '#',
            descr="Result is Valid"
        ),
        Metric(
            "SumDistance",
            "SumDistance",
            None,
            '#',
            descr="Sum Distance"
        ),
    ]

]


def parse_msi_metrics_file(f):
    """
    Nothing yet
    """

    s_name = re.search(r"(.*)\.microsat_output.json", f["fn"]).group(1)

    try:
        data = json.loads(f["f"])
    except Exception as e:
        log.debug(e)
        log.warning("Could not parse msi JSON: '{}'".format(f["fn"]))
        return

    del data["Settings"]

    return s_name, data
