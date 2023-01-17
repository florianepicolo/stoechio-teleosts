#! /usr/bin/env python3
# author : Picolo Floriane

import glob     # pour le parcours de fichiers
import re       # pour utiliser le multiple splicing
import csv      # pour utiliser l'écriture csv

def create_human_list(wholegenome_file, indice, sep):
    """
        Prend les fichiers avec les informations et récupère sous format liste les id Ensembl
    """
    list_humangene = list()
    with open(wholegenome_file) as infile:
        for row in infile.readlines()[1:]:
            list_humangene.append(row.strip().split(sep)[indice])
    list_humangene = list(set(list_humangene))
    return list_humangene

def verification_nb_genes_pathway(interactions_file):
    list_humangene = list()
    with open(interactions_file) as infile:
        for row in infile:
            list_humangene.append(row.strip().split(";")[1])
            list_humangene.append(row.strip().split(";")[2])
    list_humangene = list(set(list_humangene))
    return list_humangene

def extract_orthologs_by_species(sp_file):
    """ 
        Fichier = biomart : humangene -> ortho fish
        Prend un fichier et récupère les orthologues pr chaque gène humain
    """
    dico_ortho = dict()
    with open(sp_file) as infile:
        sp = re.split(r',|\t', next(infile).strip())[1][:-15].replace(" ", "_")
        data = infile.readlines()[1:]
        for element in data:
            row = re.split(r',|\t', element.strip())
            if len(row) == 2:
                if row[0] not in dico_ortho.keys(): 
                    dico_ortho[row[0]] = [row[1]]
                else: dico_ortho[row[0]].append(row[1])
            else: dico_ortho[row[0]] = []
    return sp, dico_ortho

def count_orthologs(dico_ortho, list_humangene):
    """
        Permet de compter le nombre de single, dupli, tripli en fonction de la demande
    """
    count_s = count_d = count_t = 0 
    for humangene, list_fishortho in dico_ortho.items():
        if humangene in list_humangene and list_fishortho != ['']:
            if len(list_fishortho) == 1: count_s += 1
            elif len(list_fishortho) == 2: count_d += 1
            elif len(list_fishortho) > 2: count_t += 1
    return count_s, count_d, count_t

def write_nbortho_by_species(filename, humanlist, dico_sp, sp_order):
    """
        Mettre le nombre d'orthologue que possède le gène humain chez l'espèce
    """
    spamwriter = csv.writer(open(filename, "w"), delimiter=';', quoting=csv.QUOTE_NONE, quotechar='"', escapechar='')
    sp_order.insert(0,"ensembl_id")
    spamwriter.writerow(sp_order)

    for humangene in humanlist:
        row = [humangene]
        for sp in sp_order[1:]: 
            if humangene in dico_sp[sp].keys():
                row.append(str(len(dico_sp[sp][humangene])))
        spamwriter.writerow(row)

def write_chi2_by_species(filename, pathwayhumanlist, wholehumanlist, dico_sp, sp_order):
    """
        Mettre le début des statistiques pour envoyer sur R par la suite
    """
    spamwriter = csv.writer(open(filename, "w"), delimiter=';', quoting=csv.QUOTE_NONE, quotechar='"', escapechar='')
    fields = ["species","nb_g_single","nb_g_dupli","nb_g_tripli", "nb_p_single","nb_p_dupli","nb_p_tripli"]
    spamwriter.writerow(fields)

    for sp in sp_order[1:]:
        gs, gd, gt = count_orthologs(dico_sp[sp], wholehumanlist)
        ps, pd, pt = count_orthologs(dico_sp[sp], pathwayhumanlist)
        spamwriter.writerow([sp, gs, gd, gt, ps, pd, pt])
        

if __name__ == "__main__":
    
    path_fish = "teleosteens/"
    file_wholegenomehuman = "g-infos-biomart.csv"
    # file_pathwayhuman = "p-infos-kegg.csv"
    # file_pathwayhuman = "p-min-interactions-ensembl.csv"
    file_pathwayhuman = "p-max-interactions-ensembl.csv"
    outfile_nbortho = "g-nbortho-teleosteens.csv"
    # outfile_chi2 = "p-min-stats-chi2-teleosteens.csv"
    outfile_chi2 = "p-max-stats-chi2-teleosteens.csv"
    ordrespecie = ["European_seabass","Large_yellow_croaker","Gilthead_seabream","Lumpfish","Stickleback","Channel_bull_blenny","Pike-perch","Ballan_wrasse","Fugu","Tetraodon","Climbing_perch","Siamese_fighting_fish","Zig-zag_eel","Barramundi_perch","Greater_amberjack","Yellowtail_amberjack","Tongue_sole","Turbot","Chinese_medaka","Japanese_medaka_HdrR","Indian_medaka","Javanese_ricefish","Turquoise_killifish","Mangrove_rivulus","Sheepshead_minnow","Mummichog","Amazon_molly","Sailfin_molly","Guppy","Platyfish","Eastern_happy","Zebra_mbuna","Makobe_Island_cichlid","Burton's_mouthbrooder","Lyretail_cichlid","Nile_tilapia","Midas_cichlid","Spiny_chromis","Clown_anemonefish","Orange_clownfish","Bicolor_damselfish","Tiger_tail_seahorse","Pinecone_soldierfish","Atlantic_cod","Northern_pike","Huchen","Chinook_salmon","Coho_salmon","Rainbow_trout","Atlantic_salmon","Brown_trout","Atlantic_herring","Denticle_herring","Mexican_tetra","Red-bellied_piranha","Electric_eel","Channel_catfish","Goldfish","Common_carp","Golden-line_barbel","Zebrafish","Paramormyrops_kingsleyae","Asian_bonytongue"]

    list_wholehumangene = create_human_list(file_wholegenomehuman, 0, "\t")
    # list_pathwayhumangene = create_human_list(file_pathwayhuman, 1, ";")
    list_pathwayhumangene = verification_nb_genes_pathway(file_pathwayhuman)

    dict_sp = dict()
    for file_fish in glob.glob(path_fish + "*.txt"):
        specie, dict_ortho = extract_orthologs_by_species(file_fish)
        dict_sp[specie] = dict_ortho
    write_nbortho_by_species(outfile_nbortho, list_wholehumangene, dict_sp, ordrespecie)
    write_chi2_by_species(outfile_chi2, list_pathwayhumangene, list_wholehumangene, dict_sp, ordrespecie)
    


