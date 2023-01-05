#! /usr/bin/env python3
# author : Picolo Floriane

import csv 
from bs4 import BeautifulSoup
import requests


def recover_missing_ensembl_id(kegg_id):
    response = requests.get("https://www.kegg.jp/entry/"+kegg_id)
    soup = BeautifulSoup(response.content, 'html.parser')
    for link in soup.find_all('a'):
        text = link.get_text('href')
        if "ENS" in text:
            return text
    return False

def get_corr_kegg_ENS():
    dico_corr = dict()
    # with open("g-ENSvsNCBI-humangenes.txt") as corrfile:
    with open("g-infos-biomart.csv") as corrfile:

        # fichier récupérer en faisant une recherche biomart des id ENS -> id NCBI (kegg = hsa:'+ id NCBI')
        for row in corrfile:
            lign = row.strip().split("\t")
            ens = lign[0]
            kegg = "hsa:"+lign[1]

            dico_corr[kegg] = ens

    return dico_corr

def create_file_interact_ensembl(dict_corr):
    with open("p-interactions-kegg.csv") as interact : 
        spamwriter = csv.writer(open("p-interactions-ensembl.csv", "w"), delimiter=';', quoting=csv.QUOTE_NONE, quotechar='"', escapechar='')
        fields = ["pathway","interact1_ensembl","interact2_ensembl"]
        spamwriter.writerow(fields)
        for row in interact:
            
            lign = row.strip().split(";")
            path = lign[0]
            i1 = lign[1]
            i2 = lign[2]

            if path != "pathway":
                if i1 in dict_corr.keys(): ens_i1 = dict_corr[i1]
                else: ens_i1 = recover_missing_ensembl_id(i1)
                if i2 in dict_corr.keys(): ens_i2 = dict_corr[i2]
                else: ens_i2 = recover_missing_ensembl_id(i2)
                l = [path, ens_i1, ens_i2]
                spamwriter.writerow(l)


def create_file_matrice_interact():
    dico_matrice = dict()

    with open("p-interactions-ensembl.csv") as interactfile :
        for row in interactfile:
            ligne = row.strip().split(";")
            path = ligne[0]
            i1 = ligne[1]
            i2 = ligne[2]

            if path not in dico_matrice.keys():
                dico_matrice[path] = [i1, i2]
                # si on veut les trier dans un ordre croissant :
                # if i1 > i2: dico_matrice[path] = [i1, i2]
                # else: dico_matrice[path] = [i2, i1]
            else: 
                dico_matrice[path].append([i1, i2])
                # if i1 > i2: dico_matrice[path].append([i1, i2])
                # else: dico_matrice[path].append([i2, i1])
    
    spamwriter = csv.writer(open("p-matrice-interactions.csv", "w"), delimiter=';', quoting=csv.QUOTE_NONE, quotechar='"', escapechar='')
    fields = dico_matrice.keys()
    spamwriter.writerow(fields)

    del(dico_matrice["pathway"])

    #     fields = ["pathway","interact1_ensembl","interact2_ensembl"]
    #     spamwriter.writerow(fields)


    for path1, list_i1 in dico_matrice.items():
        l = list()
        l.append(path1)

        for path2, list_i2 in dico_matrice.items():
            if path1 == path2: 
                l.append(len(dico_matrice[path1]))
            else: 
                nb = 0
                for couple1 in list_i1:
                    for couple2 in list_i2: 
                        if couple1 == couple2: 
                            nb += 1
                l.append(nb)

        spamwriter.writerow(l)

def create_file_matrice_genes():
    dico_matrice = dict()



if __name__ == "__main__" :
    # recover_missing_ensembl_id(kegg_id="hsa:102724428")
    # d_correspondance = get_corr_kegg_ENS()
    # create_file_interact_ensembl(d_correspondance)
    create_file_matrice_interact()


