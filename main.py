from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align.Applications import MafftCommandline
import utils, os, align
from utils import Color


def exo_a():
    request = "(SARS-CoV-2 [orgn] AND refseq [filter]) OR (Bat coronavirus RaTG13) OR (MP789 MT121216)"
    utils.request_ncbi("Nucleotide", request, "gb", "files/seq_covid.gb")

def exo_b(filepath: str, output_path: str="info_seq_covid.txt"):
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
            
            # ◦ Ecriture des informations relatif à chaque gène dans le fichier de sortie :
            for feature in feature_list:
                gene = feature.qualifiers["gene"][0]
                start = str(feature.location.start)
                end = str(feature.location.end)
                protein_id = feature.qualifiers["protein_id"][0]

                fd.write("      - Nom : " + gene + "\n")
                fd.write("        ├──Début : " + start + "\n")
                fd.write("        ├──Fin : " + end + "\n")
                fd.write("        └──Id protéine : " + protein_id + "\n")


def exo_c(gene: str, filepath: str):
    # ◦ Récupération des objets SeqRecord à partir des entrées de 'seq_covid.gb' :
    record_list = list(SeqIO.parse("files/seq_covid.gb", "genbank"))
    seqrecord_list = []
    for record in record_list:
        for feature in record.features:
            # ◦ On vérifie si le gène est bien celui recherché en paramètre de la fonction :
            if feature.type == "CDS" and feature.qualifiers["gene"][0] == gene:
                # ◦ Récupération des informations relatif au gène trouvé.
                name = feature.qualifiers["gene"][0]
                seq = Seq(feature.qualifiers["translation"][0])
                id = feature.qualifiers["protein_id"][0]
                description = feature.qualifiers["product"][0]
                # ◦ Insertion dans une liste d'un objet SeqRecord manuellement créer.
                seqrecord_list.append(SeqRecord(seq, id, name, description))
    # ◦ Ecriture du fichier FASTA grâce à la liste de SeqRecord précédemment créée.
    SeqIO.write(seqrecord_list, filepath, "fasta")


def exo_d(filepath: str, output_path: str):
    # ◦ Utilisation du wrapper d'alignement Mafft.
    command = MafftCommandline(input=filepath)
    # ◦ Éxécution de la commande.
    stdout, stderr = command()
    with open(output_path, 'w') as fd:
        # ◦ Écriture du résultat dans un fichier.
        fd.write(stdout)


def exo_e(filepath: str, output_path: str):
    # ◦ Récupération des séquences alignées
    record_list = list(SeqIO.parse(filepath, "fasta")) 
    with open(output_path, 'w') as fd:
        fd.write("position     PANGOLIN    CHAUVE-SOURIS         HOMME\n")
        for i in range(len(record_list[0].seq)):
            prot_h = record_list[0].seq[i]
            prot_c = record_list[1].seq[i]
            prot_p = record_list[2].seq[i]
            # ◦ Sélection des lignes avec des différences 
            if prot_h != prot_c or prot_h != prot_p or prot_c != prot_p:
                fd.write(f"   {i+1} " + " "*(12-len(str(i+1))) + f"{prot_h}              {prot_c}                 {prot_p}\n")


def exo_f(filepath: str):
    # ◦ Récupération des séquences et assignation de variables pour faciliter la lecture.
    record_list = list(SeqIO.parse(filepath, "fasta"))
    same_chauve, same_pangolin = 0, 0
    seq_pangolin, seq_chauve, seq_homme = record_list[0].seq, record_list[1].seq, record_list[2].seq
    for i in range(len(seq_pangolin)):
        # ◦ Incrémente de un si les deux protéines comparées sont les mêmes.
        same_chauve += 1 if seq_homme[i] == seq_chauve[i] else 0
        same_pangolin += 1 if seq_homme[i] == seq_pangolin[i] else 0
    print("  • Taux de conservation pour la chauve-souris : " + str(round((same_chauve/len(seq_chauve))*100, 1)) + "%")
    print("  • Taux de conservation pour le pangolin : " + str(round((same_pangolin/len(seq_pangolin))*100, 1)) + "%")


def analyse_gene(gene: str, name: str, type: str):
    if not name in os.listdir("files/"):
        os.mkdir(f"files/{name}")
        
    print(f"── Analyse du gène : {name}")
    
    print("Création du fichier FASTA...", end="")
    exo_c(gene, f"files/{name}/{name}.fasta")
    print(" Terminé.")
    
    print("Alignement des séquences...", end="")
    if type == "mafft":
        exo_d(f"files/{name}/{name}.fasta", f"files/{name}/aln-{name}.fasta")
    else:
        exo_k(f"files/{name}/{name}.fasta", f"files/{name}/aln-{name}.fasta")
    print(" Terminé.")
    
    print("Comparaison des séquences alignées...", end="")
    exo_e(f"files/{name}/aln-{name}.fasta", f"files/{name}/comparaison-{name}.txt")
    print(" Terminé.")
    
    print("Calcul du taux de conservation des séquences :")
    exo_f(f"files/{name}/aln-{name}.fasta")
    
    print("Conversion des gènes au format genbank")
    exo_i(gene, name)
    print(Color.GREEN + "Analyse exécutée avec succès." + Color.RESET)

def exo_i(gene: str, name: str):
    # ◦ Récupération des objets SeqRecord à partir des entrées de 'seq_covid.gb' :
    record_list = list(SeqIO.parse("files/seq_covid.gb", "genbank"))
    req_id = ""
    for record in record_list:
        for feature in record.features:
            # ◦ On vérifie si le gène est bien celui recherché en paramètre de la fonction :
            if feature.type == "CDS" and feature.qualifiers["gene"][0] == gene:
                # ◦ Récupération des id nécessaires
                req_id = req_id + feature.qualifiers["protein_id"][0] + " "
    utils.request_ncbi("Protein", req_id, "gb", f"files/{name}/out_{name}.gb")
    
def exo_j(filepath, outpath):
    record_list = list(SeqIO.parse(filepath, "fasta"))
    seqrecords = []
    #seq3=homme
    #seq2=cs
    #seq1=pangolin
    with open(outpath, "w") as fd:
        for record in [(record_list[2],record_list[1]),(record_list[2],record_list[0]),(record_list[0],record_list[1])]:
            seqa, seqb = align.pair(str(record[0].seq),str(record[1].seq))
            seqrecords.append(SeqRecord(Seq(seqa), id=record[0].id, description=record[0].description))
            seqrecords.append(SeqRecord(Seq(seqb), id=record[1].id, description=record[1].description))
    SeqIO.write(seqrecords, outpath, "fasta")

"""
Tentative de résolution pour l'étape K.
"""
def exo_k(filepath, output):
    record_list = list(SeqIO.parse(filepath, "fasta"))
    aligned_seq = []
    if len(record_list) > 1:
        first_seq = record_list[0].seq
        aligned_seq.append(first_seq)
        for i in range(1, len(record_list)):
            seq_i = record_list[i].seq
            if len(first_seq) > len(seq_i):
                _,alignement = align.pair(first_seq, seq_i)
            else:
                alignement,_ = align.pair(first_seq, seq_i)
            aligned_seq.append(Seq(alignement))
    else:
        raise "Un fichier FASTA doit contenir au minimum 2 séquences pour être aligné."
    align_record_list = []
    for i in range(len(aligned_seq)):
        og_record: SeqRecord = record_list[i]
        align_record_list.append(SeqRecord(aligned_seq[i], og_record.id, og_record.name, og_record.description))
    SeqIO.write(align_record_list, output, "fasta")

if __name__ == "__main__":
    
    #------- ◦ Décommentez le code que vous souhaitez lancer ◦ -------
    
    # exo_a()
    
    #------- ◦ Analyse automatisée ◦ -------
    
    # for gene in [{"gene":"S","name":"spike"},{"gene":"M","name":"membrane"},{"gene":"N","name":"nucleocapsid"}]:
        # print("------------------------------------")
        # print(Color.ORANGE + f"Éxécution de l'analyse du gène " + gene["name"] + " (" + gene["gene"] + ")." + Color.RESET)
        # analyse_gene(gene["gene"], gene["name"], type="mafft")
    
    #------- ◦ Etape avancées ◦ -------
    
    # exo_i("S", "spike")
    
    # exo_j("files/spike/spike.fasta", "files/spike/aln-spike-1.fasta")
    # exo_j("files/nucleocapsid/nucleocapsid.fasta", "files/nucleocapsid/aln-nucleocapsid-1.fasta")
    # exo_j("files/membrane/membrane.fasta", "files/membrane/aln-membrane-1.fasta")
    
    # exo_k("files/spike/spike.fasta", "files/spike/aln-spike-2.fasta")
    # exo_k("files/nucleocapsid/nucleocapsid.fasta", "files/nucleocapsid/aln-nucleocapsid-2.fasta")
    # exo_k("files/membrane/membrane.fasta", "files/membrane/aln-membrane-2.fasta")
    
    pass