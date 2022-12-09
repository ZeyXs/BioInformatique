from utils import Color
from main import *

def choose_etape():
    answer = input("Votre choix : " + Color.CYAN)
    print(Color.RESET)
    if not answer.isalpha():
        print("Veuillez rentrer une étape valide.")
        return choose_etape()
    return answer

def choose_gene():
    gene = input("Veuillez choisir le gène à analyser (Ex: S) : " + Color.CYAN)
    name = input(Color.RESET + "Veuillez choisir le nom du gène (Ex: spike) : " + Color.CYAN)
    print(Color.RESET)
    return gene, name

def choose_seq():
    seq1 = input("Veuillez choisir une première séquence à aligner : " + Color.CYAN)
    print(Color.RESET)
    if not seq1.isalpha():
        print(Color.RED + "Veuillez rentrer une séquence valide." + Color.RESET)
        return choose_seq()
    seq2 = input("Veuillez choisir une deuxième séquence à aligner : " + Color.CYAN)
    print(Color.RESET)
    if not seq2.isalpha():
        print(Color.RED + "Veuillez rentrer une séquence valide." + Color.RESET)
        return choose_seq
    return seq1, seq2

def select_etape():
    print("Veuillez choisir une étape à lancer : " + Color.GREEN + "(a,b,c,d,e,f,g,i,j,analyse)" + Color.RESET)
    print("PS: L'étape 'analyse' correspond à la fusion de tous les exercices qui consiste à l'analyse d'un gène.")
    etape = choose_etape()
    match etape:
        case "a":
            print(Color.ORANGE + "Éxécution de l'étape a." + Color.RESET)
            exo_a()
        case "b":
            print(Color.ORANGE + "Éxécution de l'étape b." + Color.RESET)
            exo_b("files/seq_covid.gb")
        case "c":
            print(Color.ORANGE + "Éxécution de l'étape c." + Color.RESET)
            exo_c("files/seq_covid.gb")
        case "d":
            print(Color.ORANGE + "Éxécution de l'étape d." + Color.RESET)
            exo_d("files/spike.fasta", "files/aln-spike.fasta")
        case "e":
            print(Color.ORANGE + "Éxécution de l'étape e." + Color.RESET)
            exo_e("files/comparaison-spike.txt")
        case "f":
            print(Color.ORANGE + "Éxécution de l'étape f." + Color.RESET)
            exo_f("files/spike/aln-spike.fasta")
        case "analyse":
            gene, name = choose_gene()
            print(Color.ORANGE + f"Éxécution de l'analyse du gène {name} ({gene})." + Color.RESET)
            analyse_gene(gene, name)
        case "i":
            gene, name = choose_gene()
            print(Color.ORANGE + f"Éxécution de l'étape i pour le gène {name} ({gene})." + Color.RESET)
            exo_i(gene, name)
        case "j":
            seq1, seq2 = choose_seq()
            print(Color.ORANGE + f"Éxécution de l'alignement avec les deux séquences sélectionnées." + Color.RESET)
            exo_j(seq1, seq2)
        case _:
            print(Color.RED + "Veuillez renseigner une étape valide :" + Color.RESET)
            select_etape()
            
select_etape()