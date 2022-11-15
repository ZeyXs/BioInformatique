from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align.Applications import MafftCommandline
import utils


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
                # ◦ Insertion dans une liste d'un objet SeqRecord manuellement créer :
                seqrecord_list.append(SeqRecord(seq, id, name, description))
    # ◦ Ecriture du fichier FASTA grâce à la liste de SeqRecord précédemment créée.
    SeqIO.write(seqrecord_list, filepath, "fasta")


def exo_d(filepath: str, out_path: str):
    # ◦ Utilisation du wrapper d'alignement Mafft.
    command = MafftCommandline(input=filepath)
    stdout, stderr = command()
    with open(out_path, 'w') as fd:
        fd.write(stdout)


def exo_e(out_path: str):
    record_list = list(SeqIO.parse("files/aln-spike.fasta", "fasta"))
    with open(out_path, 'w') as fd:
        fd.write("position      HOMME      CHAUVE-SOURIS       PANGOLIN\n")
        for i in range(len(record_list[0].seq)):
            fd.write(f"   {i+1} " + " "*(12-len(str(i+1))) + f"{record_list[0].seq[i]}              {record_list[1].seq[i]}                 {record_list[2].seq[i]}\n")

# TODO
def exo_f(filepath: str):
    record_list = list(SeqIO.parse(filepath, "fasta"))
    same_chauve, same_pangolin = 0, 0
    seq_homme, seq_chauve, seq_pangolin = record_list[0].seq, record_list[1].seq, record_list[2].seq
    for i in range(len(seq_homme)):
        same_chauve += 1 if seq_homme[i] == seq_chauve[i] else 0
        same_pangolin += 1 if seq_homme[i] == seq_pangolin[i] else 0
    print("Taux de conservation pour la chauve-souris : " + str(round((same_chauve/len(seq_chauve))*100, 1)) + "%")
    print("Taux de conservation pour le pangolin : " + str(round((same_pangolin/len(seq_pangolin))*100, 1)) + "%")

if __name__ == "__main__":
    #exo_a()
    #exo_b("files/seq_covid.gb")
    #exo_c("S", "files/spike.fasta")
    #exo_d("files/spike.fasta", "files/aln-spike.fasta")
    #exo_e("files/resultatComparaison_geneS.txt")
    #exo_f("files/aln-spike.fasta")
    pass