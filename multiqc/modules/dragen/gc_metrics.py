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


NAMESPACE = "gc metrics"

class DragenGCMetrics(BaseMultiqcModule):
    def add_gc_metrics(self):
        data_by_sample = defaultdict(dict)

        for f in self.find_log_files("dragen/gc_metrics"):
            s_name, data = parse_gc_metrics_file(f)
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
        self.write_data_file(data_by_sample, "dragen_gc")

        all_metric_names = set()
        for sn, sdata in data_by_sample.items():
            for m in sdata.keys():
                all_metric_names.add(m)

        gen_stats_headers, gc_table_headers = make_headers(all_metric_names, GC_METRICS)

        for s in data_by_sample:
            for m in all_metric_names:
                if m not in data_by_sample[s].keys():
                    data_by_sample[s][m] = ""

        self.general_stats_addcols(data_by_sample, gen_stats_headers, namespace=NAMESPACE)


        self.add_section(
            name="GC Metrics",
            anchor="dragen-gc",
            description="""
            Duration of the pipeline exectution.
            """,
            plot=table.plot(data_by_sample, gc_table_headers, pconfig={"namespace": NAMESPACE}),
        )
        return data_by_sample.keys()


GC_METRICS = [
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
            "Window size",
            "Window size",
            '#',
            '#',
            descr="GC Window size"
        ),
        Metric(
            "Number of valid windows",
            "Number of valid windows",
            '#',
            '#',
            descr="GC Number of valid windows"
        ),
        Metric(
            "Average reference GC",
            "Average reference GC",
            'hid',
            '#',
            descr="Average reference GC"
        ),
        Metric(
            "Mean global coverage",
            "Mean global coverage",
            '#',
            '#',
            descr="Mean global coverage"
        ),
        Metric(
            "Normalized coverage at GCs 0-19",
            "Normalized coverage at GCs 0-19",
            '#',
            '#',
            descr="Normalized coverage at GCs 0-19"
        ),
        Metric(
            "Normalized coverage at GCs 20-39",
            "Normalized coverage at GCs 20-39",
            '#',
            '#',
            descr="Normalized coverage at GCs 20-39"
        ),
        Metric(
            "Normalized coverage at GCs 40-59",
            "Normalized coverage at GCs 40-59",
            '#',
            '#',
            descr="Normalized coverage at GCs 40-59"
        ),
        Metric(
            "Normalized coverage at GCs 60-79",
            "Normalized coverage at GCs 60-79",
            '#',
            '#',
            descr="Normalized coverage at GCs 60-79"
        ),
        Metric(
            "Normalized coverage at GCs 80-100",
            "Normalized coverage at GCs 80-100",
            '#',
            '#',
            descr="Normalized coverage at GCs 80-100"
        ),
        Metric(
            "AT Dropout",
            "AT Dropout",
            '#',
            '#',
            descr="AT Dropout"
        ),
        Metric(
            "GC Dropout",
            "GC Dropout",
            '#',
            '#',
            descr="GC Dropout"
        )
    ]

]


def parse_gc_metrics_file(f):
    """
    Example.gc_metrics.csv

    GC BIAS DETAILS,,Windows at GC 0,132988,0.005
    GC BIAS DETAILS,,Windows at GC 1,95199,0.003
    GC BIAS DETAILS,,Windows at GC 2,109792,0.004
    GC BIAS DETAILS,,Windows at GC 3,140222,0.005
    GC BIAS DETAILS,,Windows at GC 4,150722,0.005
    GC BIAS DETAILS,,Windows at GC 5,153418,0.005
    GC BIAS DETAILS,,Windows at GC 6,174673,0.006
    GC BIAS DETAILS,,Windows at GC 7,188242,0.007
    GC BIAS DETAILS,,Windows at GC 8,214901,0.008
    GC BIAS DETAILS,,Windows at GC 9,243345,0.009
    GC BIAS DETAILS,,Windows at GC 10,287961,0.010
    GC BIAS DETAILS,,Windows at GC 11,353018,0.012
    GC BIAS DETAILS,,Windows at GC 12,461843,0.016
    GC BIAS DETAILS,,Windows at GC 13,610524,0.021
    GC BIAS DETAILS,,Windows at GC 14,866615,0.030
    GC BIAS DETAILS,,Windows at GC 15,1276618,0.045
    GC BIAS DETAILS,,Windows at GC 16,1911027,0.067
    GC BIAS DETAILS,,Windows at GC 17,2870773,0.100
    GC BIAS DETAILS,,Windows at GC 18,4306582,0.151
    GC BIAS DETAILS,,Windows at GC 19,6357009,0.222
    GC BIAS DETAILS,,Windows at GC 20,9183190,0.321
    GC BIAS DETAILS,,Windows at GC 21,12917454,0.451
    GC BIAS DETAILS,,Windows at GC 22,17640188,0.617
    GC BIAS DETAILS,,Windows at GC 23,23418465,0.818
    GC BIAS DETAILS,,Windows at GC 24,30230969,1.057
    GC BIAS DETAILS,,Windows at GC 25,37915855,1.325
    GC BIAS DETAILS,,Windows at GC 26,46182567,1.614
    GC BIAS DETAILS,,Windows at GC 27,54815893,1.916
    GC BIAS DETAILS,,Windows at GC 28,63524441,2.220
    GC BIAS DETAILS,,Windows at GC 29,71999791,2.516
    GC BIAS DETAILS,,Windows at GC 30,80047322,2.798
    GC BIAS DETAILS,,Windows at GC 31,87470365,3.057
    GC BIAS DETAILS,,Windows at GC 32,93977185,3.284
    GC BIAS DETAILS,,Windows at GC 33,99428234,3.475
    GC BIAS DETAILS,,Windows at GC 34,103666271,3.623
    GC BIAS DETAILS,,Windows at GC 35,106872860,3.735
    GC BIAS DETAILS,,Windows at GC 36,108913722,3.806
    GC BIAS DETAILS,,Windows at GC 37,109578590,3.830
    GC BIAS DETAILS,,Windows at GC 38,108572517,3.794
    GC BIAS DETAILS,,Windows at GC 39,106196660,3.711
    GC BIAS DETAILS,,Windows at GC 40,102845632,3.594
    GC BIAS DETAILS,,Windows at GC 41,98812829,3.453
    GC BIAS DETAILS,,Windows at GC 42,94324388,3.297
    GC BIAS DETAILS,,Windows at GC 43,89840130,3.140
    GC BIAS DETAILS,,Windows at GC 44,85792103,2.998
    GC BIAS DETAILS,,Windows at GC 45,82343539,2.878
    GC BIAS DETAILS,,Windows at GC 46,79443529,2.776
    GC BIAS DETAILS,,Windows at GC 47,76991702,2.691
    GC BIAS DETAILS,,Windows at GC 48,74420023,2.601
    GC BIAS DETAILS,,Windows at GC 49,71215463,2.489
    GC BIAS DETAILS,,Windows at GC 50,67123542,2.346
    GC BIAS DETAILS,,Windows at GC 51,62330117,2.178
    GC BIAS DETAILS,,Windows at GC 52,57349290,2.004
    GC BIAS DETAILS,,Windows at GC 53,52683083,1.841
    GC BIAS DETAILS,,Windows at GC 54,48523395,1.696
    GC BIAS DETAILS,,Windows at GC 55,44858895,1.568
    GC BIAS DETAILS,,Windows at GC 56,41330606,1.444
    GC BIAS DETAILS,,Windows at GC 57,37579098,1.313
    GC BIAS DETAILS,,Windows at GC 58,33446304,1.169
    GC BIAS DETAILS,,Windows at GC 59,28961579,1.012
    GC BIAS DETAILS,,Windows at GC 60,24460896,0.855
    GC BIAS DETAILS,,Windows at GC 61,20307479,0.710
    GC BIAS DETAILS,,Windows at GC 62,16760334,0.586
    GC BIAS DETAILS,,Windows at GC 63,13834406,0.483
    GC BIAS DETAILS,,Windows at GC 64,11421281,0.399
    GC BIAS DETAILS,,Windows at GC 65,9412260,0.329
    GC BIAS DETAILS,,Windows at GC 66,7748579,0.271
    GC BIAS DETAILS,,Windows at GC 67,6292539,0.220
    GC BIAS DETAILS,,Windows at GC 68,5031060,0.176
    GC BIAS DETAILS,,Windows at GC 69,4025404,0.141
    GC BIAS DETAILS,,Windows at GC 70,3228609,0.113
    GC BIAS DETAILS,,Windows at GC 71,2586312,0.090
    GC BIAS DETAILS,,Windows at GC 72,2092749,0.073
    GC BIAS DETAILS,,Windows at GC 73,1708384,0.060
    GC BIAS DETAILS,,Windows at GC 74,1409191,0.049
    GC BIAS DETAILS,,Windows at GC 75,1187753,0.042
    GC BIAS DETAILS,,Windows at GC 76,1025655,0.036
    GC BIAS DETAILS,,Windows at GC 77,902028,0.032
    GC BIAS DETAILS,,Windows at GC 78,802356,0.028
    GC BIAS DETAILS,,Windows at GC 79,700528,0.024
    GC BIAS DETAILS,,Windows at GC 80,597861,0.021
    GC BIAS DETAILS,,Windows at GC 81,502767,0.018
    GC BIAS DETAILS,,Windows at GC 82,422249,0.015
    GC BIAS DETAILS,,Windows at GC 83,347346,0.012
    GC BIAS DETAILS,,Windows at GC 84,269298,0.009
    GC BIAS DETAILS,,Windows at GC 85,194982,0.007
    GC BIAS DETAILS,,Windows at GC 86,129590,0.005
    GC BIAS DETAILS,,Windows at GC 87,83957,0.003
    GC BIAS DETAILS,,Windows at GC 88,57883,0.002
    GC BIAS DETAILS,,Windows at GC 89,37514,0.001
    GC BIAS DETAILS,,Windows at GC 90,24589,0.001
    GC BIAS DETAILS,,Windows at GC 91,15593,0.001
    GC BIAS DETAILS,,Windows at GC 92,9212,0.000
    GC BIAS DETAILS,,Windows at GC 93,5782,0.000
    GC BIAS DETAILS,,Windows at GC 94,3352,0.000
    GC BIAS DETAILS,,Windows at GC 95,1653,0.000
    GC BIAS DETAILS,,Windows at GC 96,986,0.000
    GC BIAS DETAILS,,Windows at GC 97,658,0.000
    GC BIAS DETAILS,,Windows at GC 98,328,0.000
    GC BIAS DETAILS,,Windows at GC 99,133,0.000
    GC BIAS DETAILS,,Windows at GC 100,339,0.000
    GC BIAS DETAILS,,Normalized coverage at GC 0,0.0093
    GC BIAS DETAILS,,Normalized coverage at GC 1,0.0118
    GC BIAS DETAILS,,Normalized coverage at GC 2,0.0163
    GC BIAS DETAILS,,Normalized coverage at GC 3,0.0390
    GC BIAS DETAILS,,Normalized coverage at GC 4,0.0498
    GC BIAS DETAILS,,Normalized coverage at GC 5,0.0531
    GC BIAS DETAILS,,Normalized coverage at GC 6,0.0777
    GC BIAS DETAILS,,Normalized coverage at GC 7,0.0744
    GC BIAS DETAILS,,Normalized coverage at GC 8,0.0724
    GC BIAS DETAILS,,Normalized coverage at GC 9,0.0713
    GC BIAS DETAILS,,Normalized coverage at GC 10,0.0832
    GC BIAS DETAILS,,Normalized coverage at GC 11,0.0984
    GC BIAS DETAILS,,Normalized coverage at GC 12,0.1167
    GC BIAS DETAILS,,Normalized coverage at GC 13,0.1354
    GC BIAS DETAILS,,Normalized coverage at GC 14,0.1497
    GC BIAS DETAILS,,Normalized coverage at GC 15,0.1720
    GC BIAS DETAILS,,Normalized coverage at GC 16,0.1926
    GC BIAS DETAILS,,Normalized coverage at GC 17,0.2130
    GC BIAS DETAILS,,Normalized coverage at GC 18,0.2391
    GC BIAS DETAILS,,Normalized coverage at GC 19,0.2565
    GC BIAS DETAILS,,Normalized coverage at GC 20,0.2732
    GC BIAS DETAILS,,Normalized coverage at GC 21,0.2851
    GC BIAS DETAILS,,Normalized coverage at GC 22,0.2967
    GC BIAS DETAILS,,Normalized coverage at GC 23,0.3089
    GC BIAS DETAILS,,Normalized coverage at GC 24,0.3192
    GC BIAS DETAILS,,Normalized coverage at GC 25,0.3297
    GC BIAS DETAILS,,Normalized coverage at GC 26,0.3410
    GC BIAS DETAILS,,Normalized coverage at GC 27,0.3518
    GC BIAS DETAILS,,Normalized coverage at GC 28,0.3638
    GC BIAS DETAILS,,Normalized coverage at GC 29,0.3769
    GC BIAS DETAILS,,Normalized coverage at GC 30,0.3893
    GC BIAS DETAILS,,Normalized coverage at GC 31,0.4057
    GC BIAS DETAILS,,Normalized coverage at GC 32,0.4227
    GC BIAS DETAILS,,Normalized coverage at GC 33,0.4444
    GC BIAS DETAILS,,Normalized coverage at GC 34,0.4671
    GC BIAS DETAILS,,Normalized coverage at GC 35,0.4920
    GC BIAS DETAILS,,Normalized coverage at GC 36,0.5205
    GC BIAS DETAILS,,Normalized coverage at GC 37,0.5540
    GC BIAS DETAILS,,Normalized coverage at GC 38,0.5915
    GC BIAS DETAILS,,Normalized coverage at GC 39,0.6334
    GC BIAS DETAILS,,Normalized coverage at GC 40,0.6775
    GC BIAS DETAILS,,Normalized coverage at GC 41,0.7233
    GC BIAS DETAILS,,Normalized coverage at GC 42,0.7713
    GC BIAS DETAILS,,Normalized coverage at GC 43,0.8219
    GC BIAS DETAILS,,Normalized coverage at GC 44,0.8714
    GC BIAS DETAILS,,Normalized coverage at GC 45,0.9182
    GC BIAS DETAILS,,Normalized coverage at GC 46,0.9590
    GC BIAS DETAILS,,Normalized coverage at GC 47,1.0011
    GC BIAS DETAILS,,Normalized coverage at GC 48,1.0462
    GC BIAS DETAILS,,Normalized coverage at GC 49,1.1023
    GC BIAS DETAILS,,Normalized coverage at GC 50,1.1697
    GC BIAS DETAILS,,Normalized coverage at GC 51,1.2577
    GC BIAS DETAILS,,Normalized coverage at GC 52,1.3637
    GC BIAS DETAILS,,Normalized coverage at GC 53,1.4839
    GC BIAS DETAILS,,Normalized coverage at GC 54,1.6123
    GC BIAS DETAILS,,Normalized coverage at GC 55,1.7528
    GC BIAS DETAILS,,Normalized coverage at GC 56,1.9135
    GC BIAS DETAILS,,Normalized coverage at GC 57,2.1115
    GC BIAS DETAILS,,Normalized coverage at GC 58,2.3651
    GC BIAS DETAILS,,Normalized coverage at GC 59,2.7088
    GC BIAS DETAILS,,Normalized coverage at GC 60,3.1461
    GC BIAS DETAILS,,Normalized coverage at GC 61,3.6758
    GC BIAS DETAILS,,Normalized coverage at GC 62,4.2709
    GC BIAS DETAILS,,Normalized coverage at GC 63,4.9099
    GC BIAS DETAILS,,Normalized coverage at GC 64,5.5641
    GC BIAS DETAILS,,Normalized coverage at GC 65,6.1859
    GC BIAS DETAILS,,Normalized coverage at GC 66,6.7744
    GC BIAS DETAILS,,Normalized coverage at GC 67,7.4192
    GC BIAS DETAILS,,Normalized coverage at GC 68,7.9938
    GC BIAS DETAILS,,Normalized coverage at GC 69,8.3632
    GC BIAS DETAILS,,Normalized coverage at GC 70,8.5408
    GC BIAS DETAILS,,Normalized coverage at GC 71,8.4447
    GC BIAS DETAILS,,Normalized coverage at GC 72,8.1165
    GC BIAS DETAILS,,Normalized coverage at GC 73,7.5178
    GC BIAS DETAILS,,Normalized coverage at GC 74,6.7449
    GC BIAS DETAILS,,Normalized coverage at GC 75,5.8159
    GC BIAS DETAILS,,Normalized coverage at GC 76,4.7610
    GC BIAS DETAILS,,Normalized coverage at GC 77,3.7338
    GC BIAS DETAILS,,Normalized coverage at GC 78,2.9013
    GC BIAS DETAILS,,Normalized coverage at GC 79,2.2334
    GC BIAS DETAILS,,Normalized coverage at GC 80,1.7270
    GC BIAS DETAILS,,Normalized coverage at GC 81,1.3309
    GC BIAS DETAILS,,Normalized coverage at GC 82,1.0265
    GC BIAS DETAILS,,Normalized coverage at GC 83,0.7808
    GC BIAS DETAILS,,Normalized coverage at GC 84,0.6333
    GC BIAS DETAILS,,Normalized coverage at GC 85,0.5569
    GC BIAS DETAILS,,Normalized coverage at GC 86,0.5335
    GC BIAS DETAILS,,Normalized coverage at GC 87,0.5057
    GC BIAS DETAILS,,Normalized coverage at GC 88,0.4511
    GC BIAS DETAILS,,Normalized coverage at GC 89,0.4063
    GC BIAS DETAILS,,Normalized coverage at GC 90,0.3451
    GC BIAS DETAILS,,Normalized coverage at GC 91,0.5081
    GC BIAS DETAILS,,Normalized coverage at GC 92,0.3287
    GC BIAS DETAILS,,Normalized coverage at GC 93,0.4121
    GC BIAS DETAILS,,Normalized coverage at GC 94,0.6699
    GC BIAS DETAILS,,Normalized coverage at GC 95,9.2463
    GC BIAS DETAILS,,Normalized coverage at GC 96,5.9908
    GC BIAS DETAILS,,Normalized coverage at GC 97,10.6073
    GC BIAS DETAILS,,Normalized coverage at GC 98,33.9120
    GC BIAS DETAILS,,Normalized coverage at GC 99,27.5943
    GC BIAS DETAILS,,Normalized coverage at GC 100,0.0759
    GC METRICS SUMMARY,,Window size,100
    GC METRICS SUMMARY,,Number of valid windows,2861317133
    GC METRICS SUMMARY,,Number of discarded windows,149384759
    GC METRICS SUMMARY,,Average reference GC,40.90
    GC METRICS SUMMARY,,Mean global coverage,9.33
    GC METRICS SUMMARY,,Normalized coverage at GCs 0-19,0.20
    GC METRICS SUMMARY,,Normalized coverage at GCs 20-39,0.46
    GC METRICS SUMMARY,,Normalized coverage at GCs 40-59,1.15
    GC METRICS SUMMARY,,Normalized coverage at GCs 60-79,5.14
    GC METRICS SUMMARY,,Normalized coverage at GCs 80-100,1.07
    GC METRICS SUMMARY,,AT Dropout,30.82
    GC METRICS SUMMARY,,GC Dropout,0.02
    """

    s_name = re.search(r"(.*)\.gc_metrics.csv", f["fn"]).group(1)

    bias_data = dict()
    summary_data = dict()


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

        if analysis == "GC METRICS SUMMARY":
            summary_data[metric] = value

    return s_name, summary_data

