from Bio import SeqRecord, SeqIO
from Bio import Entrez

Entrez.email='x@yahoo.us'
# Définition d'un mail obligatoire pour utiliser Entrez
req="(SARS-CoV-2 [orgn] AND refseq [filter]) OR (Bat coronavirus RaTG13) OR (MP789 MT121216)"
#Définition de la reqête pour gagner de la place sur l'écran.
ficXML=Entrez.esearch(db="nucleotide", term= req)
#recherche des éléments sur NCBI
dic=Entrez.read(ficXML)
idlist=dic["IdList"]
#Lecture du fichier et extraction des id
ficXML.close()
#fermeture de la recherche
ficSeq=Entrez.efetch(db='nucleotide',id=idlist, rettype='gb')
#Récupération des fichiers voulu sur NCBI 
Seqtemp=SeqIO.parse(ficSeq,'gb')
ficSeq.close()
#fermeture du efetch
SeqIO.write(Seqtemp,"seq_covid.gb",'genbank')
#Lecture des résultats et création du fichier ganbank voulu 