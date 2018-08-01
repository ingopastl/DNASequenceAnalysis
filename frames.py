class CodonFrame:
    def __init__(self, codon, percentage, property):
        self.codon = codon
        self.percentage = percentage
        self.property = property

    def __str__(self):
        return str(self.codon) + " " + str(round(self.percentage, 2)) + "% " + str(self.property)