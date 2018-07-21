class MotifFrame:
    def __init__(self, motif, index):
        self.motif = motif
        self.initial_indexes = [index]


class CodonFrame:
    def __init__(self, codon, percentage, property):
        self.codon = codon
        self.percentage = percentage
        self.property = property