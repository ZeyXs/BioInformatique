from Bio import SeqIO
from Bio import Entrez

def gathering_data():
    Entrez.email = 'bio-informatique@umontpellier.fr'
    # Définition d'un mail obligatoire pour utiliser Entrez
    req = "(SARS-CoV-2 [orgn] AND refseq [filter]) OR (Bat coronavirus RaTG13) OR (MP789 MT121216)"
    # Définition de la requête dans une variable
    fichier_xml = Entrez.esearch(db="Nucleotide", term=req)
    # Recherche des éléments sur le portail web du NCBI
    dic = Entrez.read(fichier_xml)
    id_list = dic["IdList"]
    # Lecture du fichier et extraction des id
    fichier_xml.close()
    # Fermeture de la recherche
    fichier_record = Entrez.efetch(db="Nucleotide", id=id_list, rettype="gb")
    # Récupération des fichiers voulu sur NCBI
    seq_temp = SeqIO.parse(fichier_record, "genbank")
    fichier_record.close()
    # Fermeture du efetch
    SeqIO.write(seq_temp, "seq_covid.gb", "genbank")
    # Lecture des résultats et création du fichier genbank voulu