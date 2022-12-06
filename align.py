import sys
sys.setrecursionlimit(200000000)

def affiche(tab):
    maxlen=len(str(tab[0][-1]))
    for i in range(len(tab)):
        ligne=""
        for j in tab[i]:
            j=str(j)
            n=len(j)
            if n==maxlen:
                ligne+=" "+j
            else:
                ligne+=" "+" "*(1+(maxlen-n)//2)+j
        print(ligne)
                
def get_empty_matrix(seq1, seq2):
    # ◦ Initialisation de la matrice d'origine
    matrix = [[0 for _ in range(len(seq1)+1)] for _ in range(len(seq2)+1)]
    # ◦ Remplissage des cas de base
    for x in range(0, len(matrix)):
        matrix[x][0] = -x
        for y in range(0, len(matrix[0])):
            matrix[0][y] = -y
    return matrix

def score(i, j):
    # ◦ Comparaison des deux éléments de la séquence 
    return 2 if i == j else -1

def remplir(matrix, seq1, seq2, x, y):
    # ◦ On ajoute à la position x y de la matrice le chemin avec le score le plus haut.
    score_diago = matrix[x-1][y-1] + score(seq1[y-1],seq2[x-1])
    score_vert = matrix[x-1][y] + score(seq1[y-1],"-")
    score_hori = matrix[x][y-1] + score("-",seq2[x-1])
    result = max( score_diago, score_vert , score_hori)
    return result

def create_matrix(seq1, seq2):
    # ◦ Création matrice
    matrix = get_empty_matrix(seq1, seq2)
    # ◦ Remplissage de la matrice
    for x in range(1, len(matrix)):
        for y in range(1, len(matrix[0])):
            matrix[x][y] = remplir(matrix, seq1, seq2, x, y)
    
    #affiche(matrix)
    return matrix

def parcours(matrix,x,y,seq1,seq2,tab):
    if x == 1 and y == 1:
        tab.append(seq1[0])
        tab.reverse()
        return tab
    elif seq1[x-1] == seq2[y-1] or x == y:
        tab.append(seq1[x-1])
        return parcours(matrix, x-1, y-1, seq1, seq2, tab)
    elif seq1[x-1] != seq2[y-1]:
        tab.append("-")
        if x > y:
            return parcours(matrix, x-1, y, seq1, seq2, tab)
        else:
            return parcours(matrix, x, y-1, seq1, seq2, tab)
    else:
        print("Euh... il y a quelque chose qui s'est mal passé :(")

def pair(seq1, seq2):
    matrix = create_matrix(seq1,seq2)
    chemin= parcours(matrix, len(matrix[0])-1, len(matrix)-1, seq1, seq2, [])
    return ''.join(chemin)


seq1 = "ACTG"
seq2 = "ACCTG"

print(pair(seq1, seq2))