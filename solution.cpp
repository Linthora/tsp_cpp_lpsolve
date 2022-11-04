// Auteur: DANVY Alix
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <lpsolve/lp_lib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

// Pour compiler : g++ DM_prog_lin.cpp -llpsolve55 -lcolamd -ldl -o MyExe

/** 
    /!\ J'applique une facteur 10 à la position les villes dans l'écriture du fichier .dot car cela s'avéré nécéssaire pour une meilleure visualisation des points et du parcours.
    Opérations sur le fichier .dot: 
        1) Pour visualiser directement le fichier : dot -Kneato -n -Tx11 <NameDot>.dot
        2) Pour transformer le fichier dot en pdf : dot -Kneato -n -Tpdf <fichier.dot> -o <nomPDF>.pdf
*/

/**
    Cette fonction résoud un problème de type TSP.
    Un problème de type TSP (Travelling Salesman Person) est un problème qui consiste étant donné un nombre de points dans un espace à n-dimension à trouver le chemin le plus court pour parcourir toutes les villes avant de revenir à la première.
    Ici, il s'agit du même problème pour un chemin de fibre optique entre différentes villes.

    @param file est un pointeur sur le buffer qui nous permet de lire le flux de caractères du fichier. 
    @param problemName est le nom du problème donné au début du fichier et lu précédement.
    @param comment est le commentaire lu précédement au début du fichier.
    @param dimension correspond au nombre de points qui seront considérées pendant le problème.
    @param edge_weight_type mot-clé qui décrit dans quel espace se trouvent les points.
    
    @requires dimension > 2.
    @requires problemName not to be empty.
*/
void solveTSP(ifstream* file, string problemName, string comment, unsigned int dimension, string edge_weight_type);

/**
    Cette fonction permet de remettre à 0 tout les facteurs des variables du problème donné.

    @param row est le tableau des variables du problème formulé avec l'API de lp_solve.
    @param nbVariable est le nombre de variables formulé dans le problème.
*/
void empty(REAL row[], int nbVariable);

/**
    Cette fonction permet de créer une représentation en graph 2D sous le format .dot de la solution donné.

    @param dimension correspond au nombre de points qui seront considérées pendant le problème.
    @param row est le tableau des variables du problème formulé avec l'API de lp_solve.
    @param problemName est le nom du problème donné au début du fichier et lu précédement.
    @param coord est un tableau de nombre flotants correspondant aux positions des points du problème.

    @requires dimension > 2.
    @requires problemName not to be empty.
*/
void makeDot2D(int dimension, REAL row[], string problemName, float* coord);

/**
    On pourra donner en argument le nom du fichier à étudier.
    Par défaut : "a280.tsp"
*/
int main (int argc, char* argv[]) {
    string filename = "a280.tsp";

    if(argc > 1)
        filename = argv[1];

    ifstream file(filename);

    if(!file.is_open()) {
        cout << "Impossible d'ouvrir le fichier donné." << endl;
        exit(EXIT_FAILURE);
    }

    // line sert de buffer pour parser le fichier. Les autres variable servent à stocker des éléments d'informations sur le problème fourni.
    string line, problemName, type, edge_weight_type;
    string comment = "";
    unsigned int dimension;

    // Ici on parcourt les première lignes du fichier pour récupérer les éléments d'informations du problème fourni.
    while(getline(file, line)) {
        istringstream iss(line);
        string tmp;
        iss >> tmp;
        if(tmp=="NAME") {
            iss >> tmp;
            iss >> problemName;
        }
        else if(tmp=="COMMENT") {
            iss >> tmp;
            while(iss >> tmp)
                comment = comment+" "+tmp;
            // cout << comment << endl;
        }
        else if(tmp=="TYPE") {
            iss >> tmp;
            iss >> type;
        }
        else if(tmp=="DIMENSION" || tmp=="DIMENSION:") {
            if(tmp == "DIMENSION:") //ce cas n'est pas censé se produire dans ce format de fichier mais il y a une erreur dans la synthax de la ligne DIMENSION du fichier fourni du DM.
                iss >> dimension;
            else {
                iss >> tmp;
                iss >> dimension;
            }
        }
        else if(tmp=="EDGE_WEIGHT_TYPE") {
            iss >> tmp;
            iss >> edge_weight_type;
            break;
        }
    }

    // on pourra ici ajouter d'autre type à résoudre à l'avenir
    if(type=="TSP")
        solveTSP(&file, problemName, comment, dimension, edge_weight_type);
    else 
        cout << "Ce type de résolution n'est pas encore implémenté." << endl;
    file.close();
    return 0;
}


void solveTSP(ifstream* file, string problemName, string comment, unsigned int dimension, string edge_weight_type) {
    int edgeDimension = 1;
    // Pour rester modulable à l'avenir.
    // On pourrait aussi aller plus loin en proposant une class pour chaque point en fonction de leur dimension et en y intégrant une méthode de mesure de distance entre les points.
    // On ne le fera pas ici par manque de connaissance du language c++.
    if(edge_weight_type=="EUC_2D")
        edgeDimension = 2;
    else {
        cout << "Ce type de résolution n'est pas encore implémenté." << endl;
        exit(EXIT_FAILURE);
    }
    
    cout << "Récupération des coordonnées et calcul des distances entres les points:" << endl;

    // Tableau des coordonnées des points fourni dans le problème.
    float coord[dimension][edgeDimension];

    // line sert à nouveau ici de buffer pour la lecture de fichier
    string line;

    // on saute les lignes du fichiers qui ne nous intéresse pas avant les coordonnées des points données.
    while(getline(* file, line))
        if(line=="NODE_COORD_SECTION") 
            break;
    
    // On parcourt la suite du fichier contenant les coordonnées des points.
    int counter = 0;

    while(getline(* file, line)) {
        if(line=="EOF" || counter >= dimension)
            break;
        istringstream iss(line);
        string tmp;
        iss >> tmp;

        for(int i=0;i<edgeDimension;++i) {
            iss >> coord[counter][i];  
        }
        ++counter;
    }

    // tableau qui représente les distances entre chaque points
    float dist[dimension][dimension];
    // On utilise le calcul de distance sur des points en 2 dimensions. On pourrait adapter à l'avenir cette partie pour des points en N dimensions.
    for(int i=0;i<dimension;++i)
        for(int j=i;j<dimension;++j) {
            dist[i][j] = sqrt( pow((coord[i][0]-coord[j][0]),2) + pow((coord[i][1]-coord[j][1]),2) );
            dist[j][i] = dist[i][j];
        }


    cout << "Done." << endl << endl <<"Début de la modélisation du problème:" << endl << endl;

    // On calcul le nombre de variable que va contenir notre modélisation à l'avance. Ce calcul pourrait être simplifié mais pour la lisibilité pour correcteur on le laissera tel quel.
    int nbVar = (dimension*dimension-1) + (dimension-1) ;// ((n*n)-n) : nombre de variable Xij; (n-1) : nombre de variable Ui

    // On commence par créer notre modèle avec l'API lp_solve
    lprec* lp = make_lp(0, nbVar);
    REAL row[nbVar+1];

    cout << "Contraintes -> Xij = {0,1} :" << endl;
    // Contraintes : xij = 0 ou 1
    for(counter=1;counter<=(dimension*dimension)-dimension;++counter)
        set_binary(lp, counter, true);

    // tableau qui nous permet de retenir la position de la variable Xij dans le tableau de variable de lp_solve
    int x[dimension][dimension];
    
    cout << "Done." << endl << endl << "Contraintes -> (sum xij = 1) où i de 1 à n et i!=j :" << endl;
    counter=1;
    // Contraintes : (sum xij = 1) où i de 1 à n et i!=j
    for(int i=0;i<dimension;++i) {
        empty(row,nbVar);
        for(int j=0;j<dimension;++j) {
            if(i!=j) {
                x[i][j]=counter;
                row[counter++]=1; // Cette ligne incrémente directement le compteur après l'opération sur le tableau de contraintes.
            }
            else
                x[i][j]=-1;
        }
        add_constraint(lp, row, EQ, 1);
    }

    cout << "Done." << endl << endl << "Contraintes -> (sum xij = 1) où j de 1 à n et j!=i :" << endl;
    // Contraintes : (sum xij = 1) où j de 1 à n et j!=i
    for(int j=0;j<dimension;++j) {
        empty(row,nbVar);
        for(int i=0;i<dimension;++i)
            if(i!=j)
                row[x[i][j]]=1;
        add_constraint(lp, row, EQ, 1);
    }
    
    cout << "Done." << endl << endl << "Contraintes -> " << endl << "\t 1) 1 <= ui <= n-1 où i de 2 à n :" << endl << "\t 2) ui - uj + n*xij<=n-1 où i et j de 2 à n et i!=j :" << endl;
    // Avec l'incrémentation précédente de counter, on a donc le point de départ des variables Ui
    int stratU = counter; 
    // Contraintes : 1 <= ui <= n-1 && ui - uj + n*xij<=n-1 où i et j de 2 à n et i!=j
    for(int i=2;i<=dimension;++i){
        // Contraintes : 1 <= ui <= n-1 où i de 2 à n
        // cout << "e1" << endl;
        empty(row,nbVar);
        
        set_int(lp,counter,true);
        row[counter] = 1;
        add_constraint(lp, row, LE, dimension-1);
        add_constraint(lp, row, GE, 1);

        
        // Contraintes : ui - uj + n*xij<=n-1 où i et j de 2 à n et i!=j
        for(int j=2;j<=dimension;++j) {
            if(i!=j) {
                empty(row,nbVar);
                row[counter] = 1; // ui
                row[stratU+j-2] = -1; // uj
                row[x[i-1][j-1]] = dimension;
                add_constraint(lp, row, LE, dimension-1);
            }
        }
        ++counter;
    }

    cout << "Done x2." << endl << endl << "Fonction objectif -> min : (sum( sum( Cij * Xij ) où j de 1 à n et j!=i) où i de 1 à n) :" << endl;

    // Fonction Objectif : min (sum( sum( Cij * Xij ) où j de 1 à n et j!=i) où i de 1 à n)
    // Cij représente la distance entre les points i et j.
    empty(row, nbVar);
    for(int i=0;i<dimension;++i) {
        for(int j=0;j<dimension;++j)
            if(i!=j)
                row[x[i][j]] = dist[i][j];
    }
    set_obj_fn(lp, row);

    cout << "Done." << endl << endl;

    // On écrit notre modélisation lp_solve dans un fichier .lp pour pouvoir l'étudier.
    // On converti notre nom en char* depuis une variable string pour éviter les warning durant la compilation
    char writeLpName[problemName.size() + 4]; // +4 = +3 (".lp") +1 ('\0')
    strcpy(writeLpName, (problemName+".lp").c_str());
    write_lp(lp, writeLpName);

    cout << "Modélisation terminée : début de recherche de résolution. . ."<< endl;

    // On lance la résolution et on affiche le résultats de façon claire.
    if (solve(lp) == 0) {
        cout << "\nSOLUTION de " << comment << ": Distance Minimum parcourue: " << get_objective(lp) << endl;
        get_variables(lp, row);
        counter=0;
        for (int i = 0; i < dimension; i++) {
            for(int j=0;j < dimension;j++)
                if(i!=j) {
                    if(row[counter]>0)
                        printf( "X%d-%d <---> %f\n",i+1,j+1 , row[counter]);
                    counter++;
                } 
        }
        // En plus de l'affichage du résultat dans le terminal, on sauvegarde la solution sous le format .dot qui pourra être visualisé apès.
        // Ici seulement dans le cas d'un  problème en 2 dimension.
        makeDot2D(dimension, row, problemName, (float*)coord);
    } 
    else
        cout << "Pas de solution trouvé au problème donné." << endl;

}

void empty(REAL row[], int nbVariable) {
  for (int i=1;i<=nbVariable;++i)
    row[i] = 0;
}

void makeDot2D(int dimension, REAL row[], string problemName, float* coord) {
    const string fileName = problemName + ".dot";
    const int fact = 10; // cette constante nous permet d'obtenir un rendu visuellement plus intérssant en écartant un peu plus les points.
    ofstream file(fileName);

    file << "graph {" << endl ;

    // Première boucle pour déclarer les points et leurs positions dans le fichier .dot
    for(int i=0;i<dimension;++i)
        file << '\t' << i+1 << "[pos=\"" << coord[i*2]*fact << ',' << coord[i*2+1]*fact << "\"];" << endl;
    file << endl ;

    // Deuxième pour déclarer les connexions entre les villes
    int counter = 0;
    for(int i=0;i<dimension;++i){
        for(int j=0;j<dimension;++j) {
            if(i!=j){
                if(row[counter]>0)
                    file << '\t' << i+1 << " -- " << j+1 << ';' << endl;
                counter++;
            }
        }
    }

    cout << "\nFichier " << fileName << " crée pour visualiser solution."<< endl;

    file << "}" << endl ;
    file.close();
}
