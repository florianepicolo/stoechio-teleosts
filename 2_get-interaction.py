#! /usr/bin/env python3
# author : Picolo Floriane

### programme permettant de récupérer les interactions 2 à 2 dans un graph 
### entrée : graphs 
### sortie : fichier graphs;interact1;interact2

import glob
from bs4 import BeautifulSoup
import csv

def test_type(t1, t2):
    if t1 != "groupe":
        if t2 != "groupe":
            return 0            # aucun groupe
        return 2                # groupe pour type2 mais pas pour type1
    else:
        if t2 != "groupe":
            return 1            # groupe pour type1 mais pas pour type2
        return 3                # groupe pour les 2 !

def test_nbgenes(gene1, gene2):
    if len(gene1) == 1: 
        if len(gene2) == 1:
            return 0            # 1 gène de chaque côté
        return 2                # 1 gène pour le gène 1, et plusieurs gènes pr le gène 2
    else:
        if len(gene2) == 1:
            return 1            # plusieurs gènes pour le gène 1, mais 1 gène pour le gène 1
        return 3                # plusieurs gènes pour les 2 ! 

def get_association(relations, gene1, gene2):
    if test_nbgenes(gene1, gene2) == 0:                 # si un gène de chaque côté
        relations.append([gene1, gene2])
        relations.append([gene2, gene1])
    elif test_nbgenes(gene1, gene2) == 1:               # si plusieurs gènes sous gène 1
        for i in range(len(gene1)):
            relations.append([gene1[i], gene2])
            relations.append([gene2, gene1[i]])
    elif test_nbgenes(gene1, gene2) == 2:               # si plusieurs gènes sous gène 2
        for j in range(len(gene2)):
            relations.append([gene1, gene2[j]])
            relations.append([gene2[j], gene1])
    elif test_nbgenes(gene1, gene2) == 3:               # si plusieurs gènes pour les 2 !
        for i in range(len(gene1)):
            for j in range(len(gene2)):
                relations.append([gene1[i], gene2[j]])
                relations.append([gene2[j], gene1[i]])
    return relations

def get_oneway(relations, gene1, gene2):
    if test_nbgenes(gene1, gene2) == 0:                 # si un gène de chaque côté
        relations.append([gene1, gene2])
    elif test_nbgenes(gene1, gene2) == 1:               # si plusieurs gènes sous gène 1
        for i in range(len(gene1)):
            relations.append([gene1[i], gene2])
    elif test_nbgenes(gene1, gene2) == 2:               # si plusieurs gènes sous gène 2
        for j in range(len(gene2)):
            relations.append([gene1, gene2[j]])
    elif test_nbgenes(gene1, gene2) == 3:               # si plusieurs gènes pour les 2 !
        for i in range(len(gene1)):
            for j in range(len(gene2)):
                relations.append([gene1[i], gene2[j]])
    return relations

def loop_relation(filename):
    with open(filename) as htmlfile:    
        soup = BeautifulSoup(htmlfile, 'html.parser')
        
        relations = list()
        namepath = soup.find("pathway").get("title").split(" signaling")[0]

        for relation in soup.find_all("relation"):
            
            entry1 = relation.get("entry1") # id de l'élément 1 de la relation
            entry2 = relation.get("entry2") # id de l'élément 2 de la relation

            name1 = soup.find("entry", id=entry1).get("name").split() # nom de l'élément 1 de la relation
            name2 = soup.find("entry", id=entry2).get("name").split() # nom de l'élément 2 de la relation

            type1 = soup.find("entry", id=entry1).get("type") # type de l'élément de la relation 1
            type2 = soup.find("entry", id=entry2).get("type") # type de l'élément de la relation 2

            try:
                if "binding" in relation.find("subtype").get("name") or "state change" in relation.find("subtype").get("name"):
                    # si le type de la relation est une assocition : il faut faire 1 -> 2 et 2 -> 1
                    if test_type(type1, type2) == 0:                                # SI PAS DE GROUPE !
                        relations = get_association(relations, name1, name2)
                    elif test_type(type1, type2) == 1:                              # SI GROUPE POUR 1 !
                        for element1 in soup.find("entry", id=entry1).find_all("component"):
                            name1 = soup.find("entry", id=element1.get("id")).get("name")
                            relations = get_association(relations, name1, name2)
                    elif test_type(type1, type2) == 2:                              # SI GROUPE POUR 2 !
                        for element2 in soup.find("entry", id=entry2).find_all("component"):
                            name2 = soup.find("entry", id=element2.get("id")).get("name")
                            relations = get_association(relations, name1, name2)
                    elif test_type(type1, type2) == 3:                              # SI GROUPE POUR LES DEUX ! 
                        for element1 in soup.find("entry", id=entry1).find_all("component"):
                            name1 = soup.find("entry", id=element1.get("id")).get("name")
                            for element2 in soup.find("entry", id=entry2).find_all("component"):
                                name2 = soup.find("entry", id=element2.get("id")).get("name")
                                relations = get_association(relations, name1, name2)
                else: # si le type de la relation n'est une assocition : sens unilatéral !
                    if test_type(type1, type2) == 0:                                # SI PAS DE GROUPE !
                        relations = get_oneway(relations, name1, name2)
                    elif test_type(type1, type2) == 1:                              # SI GROUPE POUR 1 !
                        for element1 in soup.find("entry", id=entry1).find_all("component"):
                            name1 = soup.find("entry", id=element1.get("id")).get("name")
                            relations = get_oneway(relations, name1, name2)
                    elif test_type(type1, type2) == 2:                              # SI GROUPE POUR 2 !
                        for element2 in soup.find("entry", id=entry2).find_all("component"):
                            name2 = soup.find("entry", id=element2.get("id")).get("name")
                            relations = get_oneway(relations, name1, name2)
                    elif test_type(type1, type2) == 3:                              # SI GROUPE POUR LES DEUX ! 
                        for element1 in soup.find("entry", id=entry1).find_all("component"):
                            name1 = soup.find("entry", id=element1.get("id")).get("name")
                            for element2 in soup.find("entry", id=entry2).find_all("component"):
                                name2 = soup.find("entry", id=element2.get("id")).get("name")
                                relations = get_oneway(relations, name1, name2)             
            except AttributeError or ImportError:
                pass
    
    # supprimer les doublons
    set_relations = []
    [set_relations.append(x) for x in relations if x not in set_relations]  

    return set_relations, namepath


def linear(gene):
    if type(gene) == list:
        return gene[0]
    return gene

def write_csv_file(dict_relations):
    spamwriter = csv.writer(open("p-interactions-kegg.csv", "w"), delimiter=';', quoting=csv.QUOTE_NONE, quotechar='"', escapechar='')
    fields = ["pathway","interact1_hsa","interact2_hsa"]
    spamwriter.writerow(fields)
    for pathname, relations in dict_relations.items():
        for r in relations:         # r est une liste
            if "cpd" not in str(r) and "undefined" not in str(r): # on vire les cpd
                r1 = linear(r[0])
                r2 = linear(r[1])
                lign = [pathname, r1, r2]
                spamwriter.writerow(lign)




if __name__ == "__main__":   

    dico_paths = dict()

    for f in glob.glob("paths/*"):
        print(f)

        dico_relation, pathname = loop_relation(filename=f)
        dico_paths[pathname] = dico_relation

    write_csv_file(dict_relations=dico_paths)
