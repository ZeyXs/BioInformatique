import sys
sys.setrecursionlimit(50000)

                
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

def parcours(matrix,x,y,seq1,seq2,align_seq1,align_seq2):
    if x == 1 and y == 1:
        align_seq1.append(seq1[0])
        align_seq1.reverse()
        align_seq2.append(seq2[0])
        align_seq2.reverse()
        return align_seq1, align_seq2
    elif seq1[x-1] == seq2[y-1] or x == y:
        align_seq1.append(seq1[x-1])
        align_seq2.append(seq2[y-1])
        return parcours(matrix, x-1, y-1, seq1, seq2, align_seq1, align_seq2)
    else: #seq1[x-1] != seq2[y-1]
        if x > y:
            align_seq1.append(seq1[x-1])
            align_seq2.append("-")
            return parcours(matrix, x-1, y, seq1, seq2, align_seq1, align_seq2)
        else:
            align_seq1.append("-")
            align_seq2.append(seq2[y-1])
            return parcours(matrix, x, y-1, seq1, seq2, align_seq1, align_seq2)

def pair(seq1, seq2):
    matrix = create_matrix(seq1, seq2)
    align_seq1, align_seq2 = parcours(matrix, len(matrix[0])-1, len(matrix)-1, seq1, seq2, [], [])
    return ''.join(align_seq1), ''.join(align_seq2) 