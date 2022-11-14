from Bio import SeqIO
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align.Applications import MafftCommandline
import utils

def exo_a():
    req = "(SARS-CoV-2 [orgn] AND refseq [filter]) OR (Bat coronavirus RaTG13) OR (MP789 MT121216)"
    utils.request_ncbi("Nucleotide", req, "gb", "files/seq_covid.gb")

#exo_a()

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
            for feature in feature_list:
                gene = feature.qualifiers["gene"][0]
                start = str(feature.location.start)
                end = str(feature.location.end)
                protein_id = feature.qualifiers["protein_id"][0]

                fd.write("      - Nom : " + gene + "\n")
                fd.write("        ├──Début : " + start + "\n")
                fd.write("        ├──Fin : " + end + "\n")
                fd.write("        └──Id protéine : " + protein_id + "\n")
                
#exo_b("files/seq_covid.gb")

def exo_c(nom_gene: str, filepath: str):
    record_list = list(SeqIO.parse("files/seq_covid.gb", "genbank"))
    seqrecord_list = []
    for record in record_list:
        for feature in record.features:
            if feature.type == "CDS" and feature.qualifiers["gene"][0] == nom_gene:
                name = feature.qualifiers["gene"][0]
                seq = Seq(feature.qualifiers["translation"][0])
                id = feature.qualifiers["protein_id"][0]
                description = feature.qualifiers["product"][0]
                seqrecord_list.append(SeqRecord(seq, id, name, description))

    SeqIO.write(seqrecord_list, filepath, "fasta")

#exo_c("S", "files/spike.fasta")

def exo_d(filepath: str, out_path: str):
    command = MafftCommandline(input=filepath)
    stdout, stderr = command()
    with open(out_path, 'w') as fd:
        fd.write(stdout)

#exo_d("files/spike.fasta", "files/aln-spike.fasta")

def exo_e(out_path: str):
    record_list = list(SeqIO.parse("files/aln-spike.fasta", "fasta"))
    with open(out_path, 'w') as fd:
        fd.write("position      HOMME      CHAUVE-SOURIS       PANGOLIN\n")
        for i in range(len(record_list[0].seq)):
            fd.write(f"   {i+1} " + " "*(12-len(str(i+1))) + f"{record_list[0].seq[i]}              {record_list[1].seq[i]}                 {record_list[2].seq[i]}\n")

#exo_e("files/resultatComparaison_geneS.txt")

def exo_f():
    record_list = list(SeqIO.parse("files/aln-spike.fasta", "fasta"))
    diff_chauve = 0
    diff_pangolin = 0
    seq_homme = record_list[0].seq
    seq_chauve = record_list[1].seq
    seq_pangolin = record_list[2].seq
    for i in range(len(seq_homme)):
        diff_chauve += 1 if seq_homme[i] != seq_chauve[i] else 0
        diff_pangolin += 1 if seq_homme[i] != seq_chauve[i] else 0
        