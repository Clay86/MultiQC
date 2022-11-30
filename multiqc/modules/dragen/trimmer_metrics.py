import logging
import re
from collections import defaultdict
from typing import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

log = logging.getLogger(__name__)


class DragenTrimmerMetrics(BaseMultiqcModule):
    def add_trimmer_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files("dragen/trimmer_metrics"):
            s_name, data = parse_trimmer_metrics_file(f)
            s_name = self.clean_s_name(s_name, f)

            if s_name in data_by_sample:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f, section="stats")
            data_by_sample[s_name] = data

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)

        if not data_by_sample:
            return set()

        table_data = DragenTrimmerMetrics.__get_table_data(data_by_sample)
        header = DragenTrimmerMetrics.get_header(table_data)

        self.add_section(
            name="Trimmer Metrics",
            anchor="trimmer-metrics",
            description="""
            Metrics on trimmed reads.
            """,
            plot=table.plot(table_data, header),
        )

        return data_by_sample.keys()

    @staticmethod
    def __get_table_data(data_by_sample: dict) -> dict:
        analysis = "TRIMMER STATISTICS"
        trimmer_data = {}
        for (
            sample,
            data,
        ) in data_by_sample.items():
            trimmer_data[sample] = {}
            for metric, stat in data[analysis].items():
                number, percentage = stat
                trimmer_data[sample][metric] = number
                if percentage:
                    trimmer_data[sample][metric + " (in %)"] = percentage

        return trimmer_data

    def get_header(trimmer_data):
        shown = ["Average input read length", "Total trimmed reads (in %)", "Total trimmed bases (in %)", "Total filtered reads (in %)", \
            "Average bases trimmed per read", "Average bases trimmed per trimmed read", "Remaining poly-G K-mers R1 3prime (in %)", \
            "Remaining poly-G K-mers R2 3prime (in %)", "Quality trimmed reads filtered R1 3prime (in %)", "Quality trimmed reads filtered R2 3prime (in %)", \
            "Adapter trimmed reads filtered R1 3prime (in %)", "Adapter trimmed reads filtered R2 3prime (in %)"]
        header = []

        first_sample = list(trimmer_data.keys())[0]
        for key in trimmer_data[first_sample].keys():
            if key in shown:
                hidden = False
            else:
                hidden = True
            if "%" in key:
                to_add = {"suffix": '%'}
            else:
                to_add = {"scale" : "Greens"}
            header.append((key, dict({"title": key, "description": key, "hidden": hidden}, **to_add)))

        return OrderedDict(header)


def parse_trimmer_metrics_file(f):
    """
    sample.trimmer_metrics.csv

    TRIMMER STATISTICS,,Total input reads,18531840
    TRIMMER STATISTICS,,Total input bases,2779776000
    TRIMMER STATISTICS,,Total input bases R1,1389888000
    TRIMMER STATISTICS,,Total input bases R2,1389888000
    TRIMMER STATISTICS,,Average input read length,150
    TRIMMER STATISTICS,,Total trimmed reads,0,0.00
    TRIMMER STATISTICS,,Total trimmed bases,0,0.00
    """

    s_name = re.search(r"(.*).trimmer_metrics.csv", f["fn"]).group(1)

    data = defaultdict(dict)
    for line in f["f"].splitlines():
        tokens = line.split(",")
        if len(tokens) == 4:
            analysis, _, metric, stat = tokens
            percentage = None
        elif len(tokens) == 5:
            analysis, _, metric, stat, percentage = tokens
        else:
            raise ValueError(f"Unexpected number of tokens in line {line}")

        try:
            stat = float(stat)
        except ValueError:
            pass
        data[analysis][metric] = (stat, percentage)

    return s_name, data
