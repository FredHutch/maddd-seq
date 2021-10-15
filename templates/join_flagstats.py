#!/usr/bin/env python3

from collections import defaultdict
import os

class FlagStat:

    def __init__(self):

        self.field_order = list()
        self.sum_i = defaultdict(int)
        self.sum_j = defaultdict(int)

    def add(self, fp):
        """Add data from a flagstats file."""

        # Read in each line
        print(f"Reading data from {fp}")
        for i, j, l in self.parse_file(fp):

            # If we haven't seen this type of output before
            if l not in self.field_order:

                # Add it to the list
                self.field_order.append(l)

            # Add in the values from this line
            self.sum_i[l] += i
            self.sum_j[l] += j

    def write(self, fp):
        """Write out the aggregate data."""

        print(f"Writing out data to {fp}")
        with open(fp, 'w') as handle:

            for l in self.field_order:

                handle.write(f"{self.sum_i[l]} + {self.sum_j[l]} {l}")

    def parse_file(self, fp):
        """
        Parse the lines from a file.
        Yield tuples of i, j, l;
        '68 + 0 supplementary' -> (68, 0, 'supplementary')
        '12478 + 0 properly paired (100.00% : N/A)' -> (12478, 0, 'properly paired (100.00% : N/A)')
        """

        with open(fp, 'r') as handle:
            for line in handle:
                i, sep, j, l = line.rstrip("\\n").split(" ", 3)
                assert sep == "+", f"Unexpected line format: {line}"
                
                try:
                    i = int(i)
                except Exception as e:
                    print(f"Unexpected line format: {line}")
                    raise e
                
                try:
                    j = int(j)
                except Exception as e:
                    print(f"Unexpected line format: {line}")
                    raise e

                yield i, j, l

flagstat = FlagStat()
for fp in os.listdir("."):
    if fp.endswith(".flagstats"):
        flagstat.add(fp)
flagstat.write("${specimen}.flagstats")
print("DONE")