def affiche(tab):
    for i in range(len(tab)):
        print(tab[i])

def get_matrix(seq1, seq2):
    # ◦ Initialisation de la matrice d'origine
    matrix = [[0 for _ in range(len(seq1)+1)] for _ in range(len(seq2)+1)]
    # ◦ Remplissage des cas de base
    for x in range(0, len(matrix)):
        matrix[x][0] = -x
        for y in range(0, len(matrix[0])):
            matrix[0][y] = -y
    return matrix

def score(i,j):
    # ◦ Comparaison des deux éléments de la séquence 
    return 2 if i == j else -1
    
def remplir(matrix,seq1,seq2,x,y):
    # ◦ On ajoute à la position x y de la matrice le chemin avec le score le plus haut.
    score_diago = matrix[x-1][y-1] + score(seq1[y-1],seq2[x-1])
    score_vert = matrix[x-1][y] + score(seq1[y-1],"-")
    score_hori = matrix[x][y-1] + score("-",seq2[x-1])
    result = max( score_diago, score_vert , score_hori)
    return result

def pair(seq1, seq2):
    # ◦ Création matrice
    matrix = get_matrix(seq1, seq2)
    # ◦ Remplissage de la matrice
    for x in range(1, len(matrix)):
        for y in range(1, len(matrix[0])):
            matrix[x][y] = remplir(matrix, seq1, seq2, x, y)
    # ◦ Recherche du chemin avec le meilleur score
    
    affiche(matrix)
    return matrix
    
        
def parcours(matrix,x,y):
    #print(x,y)
    # ◦ Recherche 
    if x == 1 and y == 1: #genre sur le diapo du cm page 52, le total du score il est de 24., 
        #alors que notre algo il trouve 31, ce qui est correct mais le cours dit autrement jsp pas pk
        # on a évidemment mis les mêmes séquences que le cours
        # et y'a pas juste un erreur dans le cours? nope sur et certain
        # car j'ai regardé une vidéo d'explication sur youtube et le mec dans la vidéo disait la même chose
        # https://www.youtube.comwatch/?v=b6AGUjqIPsA&ab_channel=Vivekanand-AlgorithmEveryDay
        # si t'as envie de watch à 21min15sec
        #Je me reconnecterai dès que je suis chez moi
        return matrix[2][2]
    elif x == 1:
        return matrix[x-1][y-1] + parcours(matrix,x,y-1)
    elif y == 1:
        return matrix[x-1][y-1] + parcours(matrix,x-1,y)
    else:
        return max(matrix[x-1][y-1] + parcours(matrix,x-1,y), matrix[x-1][y-1] + parcours(matrix,x,y-1), matrix[x-1][y-1] + parcours(matrix,x-1,y-1))

seq1 = "ACGGCT"
seq2 = "ACTGT"
tab=pair(seq1, seq2)
print(parcours(tab,len(seq1),len(seq2)))


#si même chose -> diagonale
#si diff -> gauche
#