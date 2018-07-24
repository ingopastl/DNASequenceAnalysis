import io
import sys
from frames import MotifFrame
from frames import CodonFrame

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


def get_codon_frame(rna):
    m = {"UUU": 0, "UUC": 0, "UUA": 0, "UUG": 0,
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

    codon_amino = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
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

    amino_properties = {"F": ("Aromatic", "HYDROPHOBIC"),
                        "L": ("Aliphatic", "HYDROPHOBIC"),
                        "S": ("tiny", "SMALL", "POLAR"),
                        "Y": ("Aromatic", "HYDROPHOBIC", "POLAR"),
                        "C": ("HYDROPHOBIC", "SMALL"),
                        "W": ("Aromatic", "HYDROPHOBIC", "POLAR"),
                        "P": ("SMALL"),
                        "H": ("+ve", "Charged", "Aromatic", "HYDROPHONIC", "POLAR"),
                        "Q": ("POLAR"),
                        "R": ("+ve", "Charged", "POLAR"),
                        "I": ("Aliphatic", "HYDROPHOBIC"),
                        "M": ("HYDROPHOBIC"),
                        "T": ("HYDROPHOBIC", "SMALL", "POLAR"),
                        "N": ("SMALL", "POLAR"),
                        "K": ("+ve", "Charged", "POLAR", "HYDROPHOBIC"),
                        "V": ("Aliphatic", "HYDROPHOBIC", "SMALL"),
                        "A": ("tiny", "SMALL", "HYDROPHOBIC"),
                        "D": ("-ve", "Charged", "POLAR", "SMALL"),
                        "E": ("-ve", "Charged", "POLAR"),
                        "G": ("tiny", "SMALL", "HYDROPHOBIC"),
                        "-": ("Stop")}

    codon = ""
    for i in range(0, len(rna)):  # Registra a ocorrência de todos os codons na string e armazena no dict 'm'
        codon += rna[i]
        if (len(codon) == 3):
            m[codon] += 1
            codon = ""

    frame = []
    total_c = (len(rna) - (len(rna) % 3))/3  # Total de codons na string
    for key in m:
        m[key] = (m[key]/total_c) * 100  # Coloca a porcentagem do codon no dict 'm'

        frame.append(CodonFrame(key, m[key], amino_properties[codon_amino[key]]))


    return frame


# Retorna os seis frames de codon requeridos no projeto.
def get_six_codon_frames(rna):
    rna_reverse = rna[:: -1]  # Inverte a string de rna

    frames = []
    frames.append(get_codon_frame(rna))
    frames.append(get_codon_frame(rna[1:]))
    frames.append(get_codon_frame(rna[2:]))
    frames.append(get_codon_frame(rna_reverse))
    frames.append(get_codon_frame(rna_reverse[1:]))
    frames.append(get_codon_frame(rna_reverse[2:]))

    return frames


# Não sei se é isso o que o projeto requere
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

    file = io.open("motifs.txt", "w")
    for object in motif_frame_list:
        file.write(str(object.motif) + str(object.initial_indexes) + "\n")

    return motif_frame_list


# Não sei se é isso o que o projeto requere
def get_motifs_by_length(rna, min_len, max_len):
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

    get_cg_percentage(rna)

    get_six_codon_frames(rna)

    # motifs = get_motifs_by_length(rna, 3, 12)


if __name__ == "__main__":
    main()