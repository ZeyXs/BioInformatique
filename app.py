from Bio import SeqIO
from Bio import Entrez

Entrez.email = 'x@yahoo.us'
# Définition d'un mail obligatoire pour utiliser Entrez
req = "(SARS-CoV-2 [orgn] AND refseq [filter]) OR (Bat coronavirus RaTG13) OR (MP789 MT121216)"
# Définition de la requête dans une variable
fic_xml = Entrez.esearch(db="Nucleotide", term=req)
# Recherche des éléments sur le portail web du NCBI
dic = Entrez.read(fic_xml)
id_list = dic["IdList"]
#Lecture du fichier et extraction des id
fic_xml.close()
#fermeture de la recherche
fic_seq = Entrez.efetch(db="Nucleotide", id=id_list, rettype="gb")
#Récupération des fichiers voulu sur NCBI 
seq_temp = SeqIO.parse(fic_seq, "gb")
fic_seq.close()
#fermeture du efetch
SeqIO.write(seq_temp, "seq_covid.gb", "genbank")
#Lecture des résultats et création du fichier ganbank voulu 