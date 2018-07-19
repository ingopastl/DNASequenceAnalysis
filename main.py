import io
import sys


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


def main():
    file = io.open("ZIKA_Vírus_Genoma Completo.txt", "r")
    dna = ""
    for line in file:
        if (line[0] != ">"):
            dna += line.rstrip("\n")

    rna = get_rna(dna)
    rna_reverse = rna[:: -1] # Inverte a string de rna
    print(rna)
    print(rna_reverse)
    print(get_cg_percentage(rna))


if __name__ == "__main__":
    main()
