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

def choose_align_type():
    print("Veuillez choisir le type d'algorithme pour l'alignement : (mafft, etape_j)")
    type = input("Votre choix : " + Color.CYAN)
    print(Color.RESET)
    if type != "mafft" and type != "etape_j":
        print(Color.RED + "Veuillez rentrer un choix valide." + Color.RESET)
        return choose_align_type()
    return type

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
            gene, name = choose_gene()
            print(Color.ORANGE + "Éxécution de l'étape d." + Color.RESET)
            exo_d(f"files/{name}/{name}.fasta", f"files/{name}/aln-{name}.fasta")
        case "e":
            gene, name = choose_gene()
            print(Color.ORANGE + "Éxécution de l'étape e." + Color.RESET)
            exo_e(f"files/{name}/{name}.fasta", f"files/{name}/comparaison-spike.txt")
        case "f":
            gene, name = choose_gene()
            print(Color.ORANGE + "Éxécution de l'étape f." + Color.RESET)
            exo_f(f"files/{name}/aln-{name}.fasta")
        case "analyse":
            analyse()
        case "g":
            analyse()
        case "i":
            gene, name = choose_gene()
            print(Color.ORANGE + f"Éxécution de l'étape i pour le gène {name} ({gene})." + Color.RESET)
            exo_i(gene, name)
        case "j":
            seq1, seq2 = choose_seq()
            print(Color.ORANGE + f"Éxécution de l'alignement avec les deux séquences sélectionnées." + Color.RESET)
            exo_j(seq1, seq2)
        case "k":
            for gene in [{"gene":"S","name":"spike"},{"gene":"M","name":"membrane"},{"gene":"N","name":"nucleocapsid"}]:
                print(Color.ORANGE + "Alignement du gène" + gene["name"] + "(" + gene["gene"] + ")" + Color.RESET)
                exo_k("files/" + gene["name"] + "/" + gene["name"] + ".fasta", "files/" + gene["name"] + "/aln-" + gene["name"] + "2.fasta")
        case _:
            print(Color.RED + "Veuillez renseigner une étape valide." + Color.RESET)
            print("")
            select_etape()
            
def analyse():
    type = choose_align_type()
    print(Color.GREEN + f"Vous avez sélectionné l'algorithme d'alignement {type}." + Color.RESET)
    for gene in [{"gene":"S","name":"spike"},{"gene":"M","name":"membrane"},{"gene":"N","name":"nucleocapsid"}]:
        print("------------------------------------")
        print(Color.ORANGE + f"Éxécution de l'analyse du gène " + gene["name"] + " (" + gene["gene"] + ")." + Color.RESET)
        analyse_gene(gene["gene"], gene["name"], type)
            
select_etape()