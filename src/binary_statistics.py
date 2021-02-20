# -*- coding: utf-8 -*-
# Module contains class for collecting so called binary statistics for diverse tasks.
# "Binary" means that it collects two kinds of results: positive and negative.

class BinaryStatistics:
    # Class for collecting binary statistics for diverse tasks.

    def __init__(self):
        self.positive_stat = 0
        self.negative_stat = 0
    # end def __init__

    def increment_positive(self):
        # Increase positive result counter.
        self.positive_stat += 1
    # end def increment_positive

    def increment_negative(self):
        # Increase negative result counter.
        self.negative_stat += 1
    # end def increment_positive

    def get_positive_percents(self):
        # Obtain percent of positive results.
        return round(self.positive_stat / self._get_total() * 100, 2)
    # end def get_positive_percents

    def get_negative_percents(self):
        # Obtain percent of negative results.
        return round(self.negative_stat / self._get_total() * 100, 2)
    # end def get_positive_percents

    def _get_total(self):
        # Obtain total number of processed cases.
        return self.positive_stat + self.negative_stat
    # end def _get_total

# end class BinaryStatistics
