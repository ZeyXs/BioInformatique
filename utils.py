from Bio import Entrez, SeqIO

def request_ncbi(database, request, format, filepath):
    Entrez.email = "bio-informatique@etu.umontpellier.fr"
    # Recherche des éléments sur le portail web du NCBI
    fichier_xml = Entrez.esearch(db=database, term=request)
    # Lecture du fichier et extraction des id
    dic = Entrez.read(fichier_xml)
    # Fermeture de la recherche
    fichier_xml.close()
    # Récupération des fichiers voulus sur NCBI
    fichier_record = Entrez.efetch(db="Nucleotide", id=dic["IdList"], rettype=format)
    # Lecture des résultats et création du fichier genbank voulu
    seq_temp = SeqIO.parse(fichier_record, format)
    SeqIO.write(seq_temp, filepath, format)
    # Fermeture du efetch
    fichier_record.close()