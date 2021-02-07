# -*- coding: utf-8 -*-

class BinaryStatistics:

    def __init__(self):
        self.positive_stat = 0
        self.negative_stat = 0
    # end def __init__

    def increment_positive(self):
        self.positive_stat += 1
    # end def increment_positive

    def increment_negative(self):
        self.negative_stat += 1
    # end def increment_positive

    def get_positive_percents(self):
        return round(self.positive_stat / self._get_total() * 100, 2)
    # end def get_positive_percents

    def get_negative_percents(self):
        return round(self.negative_stat / self._get_total() * 100, 2)
    # end def get_positive_percents

    def _get_total(self):
        return self.positive_stat + self.negative_stat
    # end def _get_total

# end class BinaryStatistics
