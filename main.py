import io
import sys
from frames import MotifFrame

# Transformar uma string de DNA em uma string de RNA
def get_rna(dna):
    rna = ""

    for i in range(0, len(dna)):
        if (dna[i] == 'A' or dna[i] == 'T' or dna[i] == 'C' or dna[i] == 'G'):
            if (dna[i] == 'T'):
                rna = rna + "U"
            else:
                rna = rna + dna[i]
        else:
            print("Caractere inválido")
            sys.exit()
    return rna


# Retorna a porcentagem de CG em uma cadeia de DNA/RNA
def get_cg_percentage(rna):
    string_len = len(rna)
    cg_count = 0
    for i in range(0, string_len):
        if (rna[i] == "G" or rna[i] == "C"):
            cg_count += 1
    return (cg_count/string_len) * 100


def get_codon_percentage(rna):
    frequency = {"UUU": 0, "UUC": 0, "UUA": 0, "UUG": 0,
                 "UCU": 0, "UCC": 0, "UCA": 0, "UCG": 0,
                 "UAU": 0, "UAC": 0, "UAA": 0, "UAG": 0,
                 "UGU": 0, "UGC": 0, "UGA": 0, "UGG": 0,
                 "CUU": 0, "CUC": 0, "CUA": 0, "CUG": 0,
                 "CCU": 0, "CCC": 0, "CCA": 0, "CCG": 0,
                 "CAU": 0, "CAC": 0, "CAA": 0, "CAG": 0,
                 "CGU": 0, "CGC": 0, "CGA": 0, "CGG": 0,
                 "AUU": 0, "AUC": 0, "AUA": 0, "AUG": 0,
                 "ACU": 0, "ACC": 0, "ACA": 0, "ACG": 0,
                 "AAU": 0, "AAC": 0, "AAA": 0, "AAG": 0,
                 "AGU": 0, "AGC": 0, "AGA": 0, "AGG": 0,
                 "GUU": 0, "GUC": 0, "GUA": 0, "GUG": 0,
                 "GCU": 0, "GCC": 0, "GCA": 0, "GCG": 0,
                 "GAU": 0, "GAC": 0, "GAA": 0, "GAG": 0,
                 "GGU": 0, "GGC": 0, "GGA": 0, "GGG": 0}

    codon = ""
    for i in range(0, len(rna)):
        codon += rna[i]
        if (len(codon) == 3):
            frequency[codon] += 1
            codon = ""

    total_c = (len(rna) - (len(rna) % 3))/3  # Total de codons na string

    percentage = {"UUU": 0, "UUC": 0, "UUA": 0, "UUG": 0,
                  "UCU": 0, "UCC": 0, "UCA": 0, "UCG": 0,
                  "UAU": 0, "UAC": 0, "UAA": 0, "UAG": 0,
                  "UGU": 0, "UGC": 0, "UGA": 0, "UGG": 0,
                  "CUU": 0, "CUC": 0, "CUA": 0, "CUG": 0,
                  "CCU": 0, "CCC": 0, "CCA": 0, "CCG": 0,
                  "CAU": 0, "CAC": 0, "CAA": 0, "CAG": 0,
                  "CGU": 0, "CGC": 0, "CGA": 0, "CGG": 0,
                  "AUU": 0, "AUC": 0, "AUA": 0, "AUG": 0,
                  "ACU": 0, "ACC": 0, "ACA": 0, "ACG": 0,
                  "AAU": 0, "AAC": 0, "AAA": 0, "AAG": 0,
                  "AGU": 0, "AGC": 0, "AGA": 0, "AGG": 0,
                  "GUU": 0, "GUC": 0, "GUA": 0, "GUG": 0,
                  "GCU": 0, "GCC": 0, "GCA": 0, "GCG": 0,
                  "GAU": 0, "GAC": 0, "GAA": 0, "GAG": 0,
                  "GGU": 0, "GGC": 0, "GGA": 0, "GGG": 0}

    for key in frequency:
        percentage[key] = (frequency[key]/total_c) * 100

    return percentage


def get_property_percentage(c_percentage):
    codon_map = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
                 "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
                 "UAU": "Y", "UAC": "Y", "UAA": "-", "UAG": "-",
                 "UGU": "C", "UGC": "C", "UGA": "-", "UGG": "W",
                 "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
                 "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                 "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
                 "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                 "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
                 "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                 "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
                 "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
                 "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
                 "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                 "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
                 "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

    properties = {"F": ("Aromatic", "HYDROPHOBIC"),
                  "L": ("Aliphatic", "HYDROPHOBIC"),
                  "S": ("tiny", "SMALL", "POLAR")
                  }


def get_motifs(rna, motifs_len):
    motif = ""
    motif_frame_list = []
    for i in range(0, len(rna) - (motifs_len - 1)):
        motif = rna[i:i + motifs_len]
        flag = 0

        for object in motif_frame_list:
            if (object.motif == motif):
                object.initial_indexes.append(i)
                flag = 1

        if (flag == 0):
            motif_frame_list.append(MotifFrame(motif, i))

    return motif_frame_list


def dont_know(rna, min_len, max_len):
    motif_list = []
    for i in range(min_len, max_len + 1):
        motif_list += get_motifs(rna, i)
    return motif_list

def main():
    file = io.open("ZIKA_Vírus_Genoma Completo.txt", "r")
    dna = ""
    for line in file:
        if (line[0] != ">"):
            dna += line.rstrip("\n")

    rna = get_rna(dna)
    rna_reverse = rna[:: -1]  # Inverte a string de rna
    print(rna)
    print(get_cg_percentage(rna))
    codon_per = get_codon_percentage(rna)
    get_property_percentage(codon_per)

    motifs = dont_know(rna, 3, 12)

    file = io.open("motifs.txt", "w")
    for object in motifs:
        file.write(str(object.motif) + str(object.initial_indexes) + "\n")


if __name__ == "__main__":
    main()