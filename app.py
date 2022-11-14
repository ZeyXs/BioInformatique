from Bio import SeqIO
from Bio import Entrez

def gathering_data():
    # Définition d'un mail obligatoire pour utiliser Entrez
    Entrez.email = 'bio-informatique@umontpellier.fr'
    # Définition de la requête dans une variable
    req = "(SARS-CoV-2 [orgn] AND refseq [filter]) OR (Bat coronavirus RaTG13) OR (MP789 MT121216)"
    # Recherche des éléments sur le portail web du NCBI
    fichier_xml = Entrez.esearch(db="Nucleotide", term=req)
    # Lecture du fichier et extraction des id
    dic = Entrez.read(fichier_xml)
    id_list = dic["IdList"]
    # Fermeture de la recherche
    fichier_xml.close()
    # Récupération des fichiers voulu sur NCBI
    fichier_record = Entrez.efetch(db="Nucleotide", id=id_list, rettype="gb")
    # Lecture des résultats et création du fichier genbank voulu
    seq_temp = SeqIO.parse(fichier_record, "genbank")
    SeqIO.write(seq_temp, "files/seq_covid.gb", "genbank")
    # Fermeture du efetch
    fichier_record.close()

gathering_data()

def list_annotations(filepath: str, output_path: str="info_seq_covid.txt"):
    with open(f"files/{output_path}", 'w') as fd:
        record_list = list(SeqIO.parse(f"{filepath}","gb"))
        for record in record_list:
            annotations = record.annotations
            # ◦ Le nom de l’organisme dont provient la donnée ainsi que la taxonomie correspondante :
            fd.write(annotations["organism"] + " (" + ', '.join(annotations["taxonomy"]) + ")\n")
            # ◦ Le numéro d’accession de la donnée GenBank :
            fd.write("   - Numero d'accession : " + annotations["accessions"][0] + "\n")
            # ◦ La date de création de la donnée GenBank :
            fd.write("   - Date de création : " + annotations["date"] + "\n")    
            
            # Récupération de toutes les features :
            feature_list = []
            for feature in record.features:
                if feature.type == "CDS":
                    feature_list.append(feature)
                    
            # ◦ Le nombre de gènes présents sur la séquence :
            fd.write("   - Nombre de gènes présents : " + str(len(feature_list)) + "\n")
            
            # ◦ Le % de GC de la séquence du génome
            
            rv_compl = record.seq.reverse_complement()
            nb_g, nb_c = rv_compl.count("G"), rv_compl.count("C")
            fd.write("   - Pourcentage de GC : " + str(round(((nb_g + nb_c)/len(rv_compl))*100, 2)) + "%\n")
            
            # ◦ La liste des noms de gènes présents, leur position de début et de fin du CDS correspondant,
            # ainsi que l’identifiant de la protéine codée :
            fd.write("   - Liste des gènes : " + "\n")
            for feature in feature_list:
                gene = feature.qualifiers["gene"][0]
                start = str(feature.location.start)
                end = str(feature.location.end)
                protein_id = feature.qualifiers["protein_id"][0]

                fd.write("      - Nom : " + gene + "\n")
                fd.write("        ├──Début : " + start + "\n")
                fd.write("        ├──Fin : " + end + "\n")
                fd.write("        └──Id protéine : " + protein_id + "\n")


                
list_annotations("files/seq_covid.gb")