/*

Discrepancy _ Gotta go fast.
© Quentin FERRE ©
2015
**************************************************************************************************************************************************************************
Etat de la version actuelle:

1D:ok
2D:ok
3D:ok

sequence:
1D:ok
2D:ok
3D:ok

**************************************************************************************************************************************************************************

Permet, à partir d'une suite de point de calculer la discrepance.
Dimensions implémentées : 1,2 et 3D.
Possibilité d'obtenir le résultat sous forme de séquence.

**************************************************************************************************************************************************************************

Fonctions importantes:
    -void quickSort(double tableau[], int debut, int fin):
        permet de trier un tableau de dimension 1 dans l'ordre croissant.

    -void import_options(int argc, char** argv):
        lit les informations passées dans la ligne de commande et les interprète.

    -bool import_data_xD():
        importe les données de dimension x et retourne true si les données ont bien été chargées.

    -double calcul_discrepancy_xD():
        retourne la discrepance de la suite, parallélisée au niveau le plus bas(dernière boucle des axes).

    -void export_data(double):
        exporte le résultat du calcul de discrepance dans le fichier de sortie renseignée dans les paramètres.

    -sequence_ ... :
        certaines de ces fonctions se retrouvent avec sequence_ suivi du nom de la fonction: ce sont des fonctions avec le même objectif mais
            adaptées pour le traitement des données en vue du calcul de la séquence.

**************************************************************************************************************************************************************************

Variables importantes:

    -Points_xD *points_tab_xD:
        tableau de points de dimension x. Contient tous les points dans l'ordre du fichier d'entrée.

    -int pts_number:
        renseigne sur le nombre de points de la suite.

    -std::string fn_input:
        nom du fichier d'entrée.

    -std::string fn_output:
        nom du fichier de sortie.

    -int dimension:
        dimension dans laquelle on souhaite que le programme évolue.

    -double *tab_x:
        tableau 1D contenant la liste des coordonnées de l'axe x triées dans l'ordre croissant.

    -double *pos_?_point_axe_?:
        suite de tableaux permettant à partir de l'index d'une des coordonnées de retrouver les autres corrdonnées du point correspondant.

    -int* tmp_ligne:
        permet sur le plan 2D <i,j> de retenir le nombre de points dans la boite [i-1,j-1].

    -int* tmp_ligne_maj:
        permet de mettre à jour tmp_ligne.

    -int **tmp_box:
        permet, pour une boite dans un repere 3d <i,j,k>, de connaitre le nombre de point dans la boite [i,j,k-1].

    -int **tmp_box_maj:
        permet de mettre à jour tmp_box.

    -bool sequence:
        permet de savoir si l'option sequence est activée.

    -int nb_threads:
        permet de savoir sur combien de threads on souhaite executer le programme.

**************************************************************************************************************************************************************************

Méthode pour l'ajout de dimensions:

    Si vous souhaitez ajouter le traitement d'une nouvelle dimension à ce programme, les conseils suivant vous seront utiles:

        -Pour ajouter une dimension N, reprendre les fonctions N-1 et faire attention aux éléments suivant:

            -Modifier la fonction import_options pour quelle puisse traiter les nouvelles informations si vous souhaitez en ajouter.

            -Ajouter un tableaux 1D contenant les coordonnées des points sur l'axe ajouté et le trier à l'aide de la fonction quickSort.

            -Ajouter une fonction import en reprenant les modèles pour pouvoir traiter un nouveau format de données.

            -Ajouter un tableau de dimension D-1, permettant de stocker le nombre de points à l'étape précédente dans le calcul de discrepance.

            -Paralléliser les affectation des nouveaux tableaux.

            -Ajouter et initialiser les tableaux pos_?_point_axe_? permettant à partir de l'index d'une des coordonnées de retrouver les autres corrdonnées du point correspondant.

            -Dans le calcul de discrepance, modifier le calcul des différents angles de la boite.

            -A la fin de la fonction de calcul de discrepance, ajouter les calculs sur les nouveaux axes pour les boites vides.

            -Ajouter dans l'option de parallélisation firstprivate les nouveaux tableaux et variables.

            -Créer une classe Points_xD en suivant les modèles.

        -Pour ajouter l'option sequence sur la dimension souhaitée, reprendre le traitement réalisé dans le main pour l'initialisation et la mise à jour des différents tableaux
            et recréer une fonction de calcul de discrepance en suivant les modèles. Il s'agit d'apporter de petites modification à la version non séquentialisée de la fonction
            pour permettre le traitement par paliers.

**************************************************************************************************************************************************************************

*/

#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <limits>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <sys/resource.h>
#include <boost/program_options.hpp>
#include <vector>
#include <omp.h>

#include "Points_1D.h"
#include "Points_2D.h"
#include "Points_3D.h"

namespace boostPO = boost::program_options;


    Points_1D *points_tab_1D ;
    Points_2D *points_tab_2D ;
    Points_3D *points_tab_3D ;
    int pts_number;
	std::string fn_input;
	std::string fn_output;
	int dimension;

	time_t td,tf;

	double *tab_x;
	double *tab_y;
	double *tab_z;

    double *pos_y_point_axe_x;
    double *pos_y_point_axe_z;

    double *pos_x_point_axe_y;
    double *pos_x_point_axe_z;

    double *pos_z_point_axe_x;
    double *pos_z_point_axe_y;

    int* tmp_ligne;
    int* tmp_ligne_maj;

    int **tmp_box;
    int **tmp_box_maj;

    bool sequence;
    int nb_threads;



double box_area_calcul(double x, double y){
	return x*y;
}

double box_area_calcul(double x, double y,double z){
	return x*y*z;
}

void echanger(double tab[], int a, int b){
//echange la position de deux points dans un tableau, utilisé pour quicksort
	double temp = tab[a];
	tab[a] = tab[b];
	tab[b] = temp;
}

void quickSort(double tableau[], int debut, int fin){
//tri croissant du tableau
	int gauche = debut-1;
	int droite = fin+1;
	const double pivot = tableau[debut];

	if(debut >= fin)
		return;
	while(1)
	{
		do droite--; while(tableau[droite] > pivot);
		do gauche++; while(tableau[gauche] < pivot);

		if(gauche < droite)
			echanger(tableau, gauche, droite);
		else break;
	}
	quickSort(tableau, debut, droite);
	quickSort(tableau, droite+1, fin);
}

void import_options(int argc, char** argv){

	boostPO::variables_map vm;
	boostPO::options_description desc("Allowed options");
	desc.add_options()
		("help,h",
			"produce help message")
		("input,i",
			boostPO::value(&fn_input)->required(),
			"input's filename")
        ("output,o",
			boostPO::value(&fn_output)->required(),
			"output's filename")
        ("dimension,d",
			boostPO::value(&dimension)->required(),
			"dimension")
        ("sequence,s",
			"sequence enable")
        ("threads,t",
        	boostPO::value(&nb_threads)->default_value(omp_get_max_threads()),
			"number of threads")
		;

	boostPO::positional_options_description p;
	p.add("input", -1);

	try
	{
		boostPO::store(
			boostPO::command_line_parser(argc, argv).
				options(desc).positional(p).run(), vm);
		boostPO::notify(vm);
	}
	catch(boost::program_options::error& e)
	{
		std::cerr << desc << "\n";
		exit(EXIT_FAILURE);
	}

	if(vm.count("help"))
	{
		std::cout << desc << "\n";
		exit(EXIT_SUCCESS);
	}

	if(vm.count("sequence"))
	{
		std::cout << "\n Sequence : ON" << "\n";
		sequence = true;
	}
	else
	{
        std::cout << "\n Sequence : OFF" << "\n";
        sequence = false;
	}

		if(nb_threads>omp_get_max_threads())
		{
            std::cerr << "Error Threads : Number of threads available : " <<omp_get_max_threads()<< "\n";
            exit(EXIT_FAILURE);
		}
		else
            std::cout << "\n Threads : " <<nb_threads<< "\n";

    std::cout << " Dimension : " <<dimension<< "\n";
}

        ///////////////////////////  1D  ///////////////////////////////////

bool import_data_1D(){

    std::string line;

	std::ifstream fichier(fn_input.c_str(), std::ios::in);
    	if(!fichier.fail())
	{
        std::cout << "  ***************************************  \n\n" ;

		std::cout << " Chargement des donnees...\n" ;
		pts_number = 0;

		while (std::getline(fichier, line)){
		if(line != "" && line !="\n" && line!=" "){
        pts_number++;
        }
        else continue;
	}

		std::cout <<" "<<pts_number<<" lignes. \n";
		if(pts_number>0){

			fichier.clear();
			fichier.seekg(0, std::ios::beg);

			double x=0;
			int i=0;

			points_tab_1D = new Points_1D[pts_number];
			tab_x=new double[pts_number];

            while(fichier >> x )
            {
                points_tab_1D[i].set_pos_x(x);

				tab_x[i]=x;

				i++;
            }
			std::cout << "\n Data Ok !\n\n" ;
			}

			fichier.close();
			return true;
    }
    else {
    std::cerr << " Can't open input !\n" ;
    return false;
    }

}

double calcul_discrepancy_1D(){

////////////////////////// INITIALISATION /////////////////////////////////
                std::cout << " Calcul discrepance ...\n\n" ;
                double discrepancy=0;
                const int const_pt_num = pts_number;

                quickSort(tab_x,0,const_pt_num-1);

////////////////////////// CALCUL DISCREPANCE /////////////////////////////////
                #pragma omp parallel for num_threads(nb_threads)
				for(int j=0;j<const_pt_num;j++)
				{
                    int memoire=j;

                    double tmp=fabs(tab_x[j]-((2*(j+1)-1)/(double)(2*const_pt_num)));

                    #pragma omp flush(discrepancy)
                    if(discrepancy<=tmp)
                    {
                        #pragma omp crititcal
                        if(discrepancy<=tmp)
                        {
                            discrepancy=tmp;
                        }
                    }
				}

                discrepancy=discrepancy+((double)1/(double)(2*pts_number));

////////////////////// FIN DU CALCUL //////////////////////////////
                return discrepancy;
}

void export_data_1D(double discrepancy){

            std::ofstream file(fn_output.c_str());
            if(!file)
            std::cerr << "Can't open ouptut" << std::endl;

            // file << pts_number << "  " << discrepancy << std::endl;
            file  << "  " << discrepancy << std::endl;
            file.close();
}

        ///////////////////////////  2D  ///////////////////////////////////

bool import_data_2D(){

    std::string line;
	std::ifstream fichier(fn_input.c_str(), std::ios::in);

    if(!fichier.fail())
	{
        std::cout << "  ***************************************  \n\n" ;
		std::cout << " Chargement des donnees...\n" ;
		pts_number = 0;

		while (std::getline(fichier, line)){
            if(line != "" && line !="\n" && line!=" "){
                pts_number++;
            }
        else continue;
        }

		std::cout <<" "<<pts_number<<" lignes. \n";
		if(pts_number>0){

			fichier.clear();
			fichier.seekg(0, std::ios::beg);

			double x=0;
			double y=0;
			int i=0;

			points_tab_2D = new Points_2D[pts_number];
			tab_x=new double[pts_number];
			tab_y=new double[pts_number];
			pos_y_point_axe_x=new double[pts_number];
			pos_x_point_axe_y=new double[pts_number];

			tmp_ligne=new int[pts_number];
			tmp_ligne_maj=new int[pts_number];

            while(fichier >> x >> y)
            {
                points_tab_2D[i].set_pos_x(x);
				points_tab_2D[i].set_pos_y(y);

				tab_x[i]=x;
				tab_y[i]=y;
                
                

				i++;
            }

			std::cout << "\n Data Ok !\n\n" ;

        }

			fichier.close();
			return true;
    }
    else {
        std::cerr << " Can't open input !\n" ;
        return false;
    }

}

double calcul_discrepancy_2D(){

			std::cout << " Calcul discrepance ...\n\n" ;
            double discrepancy=0;

			quickSort(tab_x,0,pts_number-1);
			quickSort(tab_y,0,pts_number-1);

////////////////////////// INITIALISATION TAB POS///////////////////////////////

             for(int i=0;i<pts_number;i++)
            {
                #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<pts_number;j++)
                {
                    if(std::fabs(tab_x[i]-points_tab_2D[j].get_pos_x())<std::numeric_limits<double>::epsilon())
                        pos_y_point_axe_x[i]=points_tab_2D[j].get_pos_y();

                    if(std::fabs(tab_y[i]-points_tab_2D[j].get_pos_y())<std::numeric_limits<double>::epsilon())
                        pos_x_point_axe_y[i]=points_tab_2D[j].get_pos_x();
                }
            }

///////////////////////////////////////////////////////////////////////////

            #pragma omp parallel for num_threads(nb_threads)
			for(int i=0;i<pts_number;i++)
			{
                tmp_ligne[i]=0;
                tmp_ligne_maj[i]=0;
			}

            const int const_pt_num = pts_number;

////////////////////////// DISCREPANCE //////////////////////////////////////////////////

			for(int i=0;i<pts_number;i++)
			{
                #pragma omp parallel for num_threads(nb_threads)
                for(int k=0;k<const_pt_num;k++)
                {
                    tmp_ligne[k]=tmp_ligne_maj[k];
                }

                #pragma omp parallel for firstprivate(tmp_ligne,pos_y_point_axe_x,pos_x_point_axe_y,tab_x,tab_y) num_threads(nb_threads)
				for(int j=0;j<const_pt_num;j++)
				{
                    int memoire=0;

                    if(j>0)
                            memoire=tmp_ligne[j-1];

                    if(tab_y[j]>pos_y_point_axe_x[i])
                        {
                            memoire++;
                        }

                    if(tab_x[i]>pos_x_point_axe_y[j])
                        {
                            memoire++;
                        }

                    if(std::fabs(tab_x[i]-pos_x_point_axe_y[j])<std::numeric_limits<double>::epsilon() && std::fabs(tab_y[j]-pos_y_point_axe_x[i])<std::numeric_limits<double>::epsilon())
                        {
                            memoire++;
                        }

                    double tmp=0;
                    tmp=fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(tab_x[i],tab_y[j]));

                    if(i<const_pt_num-1 && j<const_pt_num-1)
                    {
                        tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(tab_x[i+1],tab_y[j+1])));
                        tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(tab_x[i+1],tab_y[j])));
                        tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(tab_x[i],tab_y[j+1])));
                    }

                    else if (i==const_pt_num-1 && j==const_pt_num-1)
                    {
                        tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(1,1)));
                    }

                    else if(i==const_pt_num-1)
                        tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(1,tab_y[j])));

                    else if (j==const_pt_num-1)
                        tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(tab_x[i],1)));

                    else
                        std::cout<<"impossible\n";

                    tmp_ligne_maj[j]=memoire;

                    #pragma omp flush(discrepancy)
                    if(discrepancy<=tmp)
                    {
                        #pragma omp crititcal
                        if(discrepancy<=tmp)
                        {
                            discrepancy=tmp;
                        }
                    }
				}
            }


/////////// BOUCLES SUPPLEMENTAIRES POUR TESTER LES BOITES VIDES SUR LES AXES X ET Y///////////////////////////////

            #pragma omp parallel for num_threads(nb_threads)
            for(int i=0;i<pts_number;i++)
            {
                double tmp=box_area_calcul(tab_x[i],tab_y[0]);
                #pragma omp flush(discrepancy)
                if(discrepancy<=tmp)
                {
                    #pragma omp crititcal
                    if(discrepancy<=tmp)
                    {
                        discrepancy=tmp;
                    }
                }
            }

            #pragma omp parallel for num_threads(nb_threads)
            for(int i=1;i<pts_number;i++)
            {
                double tmp=box_area_calcul(tab_x[0],tab_y[i]);
                #pragma omp flush(discrepancy)
                if(discrepancy<=tmp)
                {
                    #pragma omp crititcal
                    if(discrepancy<=tmp)
                    {
                        discrepancy=tmp;
                    }
                }
            }

////////////////////// FIN DU CALCUL //////////////////////////////

        return discrepancy;
}

void export_data_2D(double discrepancy){

            std::ofstream file(fn_output.c_str());
            if(!file)
            std::cerr << "Can't open ouptut" << std::endl;

            file << pts_number << "  " << discrepancy << std::endl;
            file.close();
}

        //////////////////////////  3D  ////////////////////////////////////

bool import_data_3D(){

        std::string line;

        std::ifstream fichier(fn_input.c_str(), std::ios::in);
    	if(!fichier.fail())
        {
            std::cout << "  ***************************************  \n\n" ;
            std::cout << " Chargement des donnees...\n" ;
            pts_number = 0;

            while (std::getline(fichier, line)){
                if(line != "" && line !="\n" && line!=" "){
                    pts_number++;
                }
            else continue;
            }

            std::cout <<" "<<pts_number<<" lignes. \n";
            if(pts_number>0){

                fichier.clear();
                fichier.seekg(0, std::ios::beg);

                double x=0;
                double y=0;
                double z=0;
                int i=0;

                points_tab_3D = new Points_3D[pts_number];
                tab_x=new double[pts_number];
                tab_y=new double[pts_number];
                tab_z=new double[pts_number];

                pos_y_point_axe_x=new double[pts_number];
                pos_y_point_axe_z=new double[pts_number];

                pos_x_point_axe_y=new double[pts_number];
                pos_x_point_axe_z=new double[pts_number];

                pos_z_point_axe_x=new double[pts_number];
                pos_z_point_axe_y=new double[pts_number];

                tmp_ligne=new int[pts_number];
                tmp_ligne_maj=new int[pts_number];

                tmp_box=new int*[pts_number];
                tmp_box_maj=new int*[pts_number];

                for(int i=0;i<pts_number;i++)
                {
                    tmp_box[i]=new int [pts_number];
                    tmp_box_maj[i]=new int [pts_number];

                }

                while(fichier >> x >> y >> z)
                {
                    points_tab_3D[ i ].set_pos_x( x );
                    points_tab_3D[ i ].set_pos_y( y );
                    points_tab_3D[ i ].set_pos_z( z );

                    tab_x[ i ] = x;
                    tab_y[ i ] = y;
                    tab_z[ i ] = z;

                    i++;
                }

                std::cout << "\n Data Ok !\n\n" ;

            }

			fichier.close();
			return true;
    }
    else {
        std::cerr << " Can't open input !\n" ;
        return false;
    }

}

double calcul_discrepancy_3D(){

			std::cout << " Calcul discrepance ...\n\n" ;
            double discrepancytmp=0;

			quickSort(tab_x,0,pts_number-1);
			quickSort(tab_y,0,pts_number-1);
			quickSort(tab_z,0,pts_number-1);

////////////////////////// INITIALISATION TAB POS///////////////////////////////

            bool x,y,z;
            double* tmp_x = new double[pts_number];
            double* tmp_y = new double[pts_number];
            double* tmp_z = new double[pts_number];

            for(int i=0;i<pts_number;i++)
            {
                tmp_x[i]=0;
                tmp_y[i]=0;
                tmp_z[i]=0;
            }


           for(int i=0;i<pts_number;i++)
            {
                x=true;
                y=true;
                z=true;
                #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<pts_number;j++)
                {
                    if(x==true){
                        if(std::fabs(tab_x[i]-points_tab_3D[j].get_pos_x())<std::numeric_limits<double>::epsilon() && tmp_x[j]==0)
                        {
                            pos_y_point_axe_x[i]=points_tab_3D[j].get_pos_y();
                            pos_z_point_axe_x[i]=points_tab_3D[j].get_pos_z();
                            x=false;
                            tmp_x[j]=5;
                        }
                    }
                    if(y==true){
                        if(std::fabs(tab_y[i]-points_tab_3D[j].get_pos_y())<std::numeric_limits<double>::epsilon() && tmp_y[j]==0)
                        {
                            pos_x_point_axe_y[i]=points_tab_3D[j].get_pos_x();
                            pos_z_point_axe_y[i]=points_tab_3D[j].get_pos_z();
                            y=false;
                            tmp_y[j]=5;
                        }
                    }
                    if(z==true){
                        if(std::fabs(tab_z[i]-points_tab_3D[j].get_pos_z())<std::numeric_limits<double>::epsilon() && tmp_z[j]==0)
                        {
                            pos_x_point_axe_z[i]=points_tab_3D[j].get_pos_x();
                            pos_y_point_axe_z[i]=points_tab_3D[j].get_pos_y();
                            z=false;
                            tmp_z[j]=5;
                        }
                    }
                }
            }

//////////////////////////////////////////////////////////////////////////

            #pragma omp parallel for num_threads(nb_threads)
			for(int i=0;i<pts_number;i++)
			{
                tmp_ligne[i]=0;
                tmp_ligne_maj[i]=0;
			}

/////////////////////////////////////////////////////////////////////////////////////////

            for(int i=0;i<pts_number;i++)
            {
              #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<pts_number;j++)
                {
                    tmp_box[i][j]=0;
                    tmp_box_maj[i][j]=0;
                }
            }

            const int test = pts_number;
////////////////////////// DISCREPANCE //////////////////////////////////////////////////

    for(int k=0; k<pts_number;k++){

            for(int o=0; o<pts_number;o++){
                #pragma omp parallel for
                for(int p=0; p<pts_number;p++){
                    tmp_box[o][p]=tmp_box_maj[o][p];

                }
            }
                #pragma omp parallel for num_threads(nb_threads)
                for(int l=0;l<pts_number;l++)
                {
                    tmp_ligne_maj[l]=0;
                }

			for(int i=0;i<pts_number;i++)
			{
                #pragma omp parallel for num_threads(nb_threads)
                for(int l=0;l<pts_number;l++)
                {
                    tmp_ligne[l]=tmp_ligne_maj[l];
                }


                #pragma omp parallel for firstprivate(tmp_ligne,tmp_box,pos_z_point_axe_x,pos_z_point_axe_y,pos_y_point_axe_x,pos_y_point_axe_z,pos_x_point_axe_y,pos_x_point_axe_z,tab_x,tab_y,tab_z) num_threads(nb_threads)
				for(int j=0;j<test;j++)
				{

                   int memoire=0;
                   int etage_inf=0;

                    if(j>0 )
                    {
                        memoire=tmp_ligne[j-1];
                    }

                    etage_inf=tmp_box[i][j];

                    if(tab_y[j]>pos_y_point_axe_x[i] && std::fabs(tab_z[k]-pos_z_point_axe_x[i])<std::numeric_limits<double>::epsilon())//////////////////////////////
                        {
                            memoire++;
                        }

                    if(tab_x[i]>pos_x_point_axe_y[j] && std::fabs(tab_z[k]-pos_z_point_axe_y[j])<std::numeric_limits<double>::epsilon())//////////////////////////////
                        {
                            memoire++;
                        }
                //si il y a un point dans la boîte
                if(
                    std::fabs(tab_x[i]-pos_x_point_axe_y[j])<std::numeric_limits<double>::epsilon() &&
                    std::fabs(tab_x[i]-pos_x_point_axe_z[k])<std::numeric_limits<double>::epsilon() &&

                    std::fabs(tab_y[j]-pos_y_point_axe_x[i])<std::numeric_limits<double>::epsilon() &&
                    std::fabs(tab_y[j]-pos_y_point_axe_z[k])<std::numeric_limits<double>::epsilon() &&

                    std::fabs(tab_z[k]-pos_z_point_axe_x[i])<std::numeric_limits<double>::epsilon() &&
                    std::fabs(tab_z[k]-pos_z_point_axe_y[j])<std::numeric_limits<double>::epsilon()
                )
                {
                    memoire++;
                }


                    double tmp=0;

                    tmp=fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],tab_y[j],tab_z[k]));
                //on check les autres angles de la boites
                   if(i<test-1 && j<test-1 && k<test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i+1],tab_y[j],tab_z[k])));
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],tab_y[j+1],tab_z[k])));
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i+1],tab_y[j+1],tab_z[k])));
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],tab_y[j],tab_z[k+1])));
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i+1],tab_y[j],tab_z[k+1])));
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],tab_y[j+1],tab_z[k+1])));
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i+1],tab_y[j+1],tab_z[k+1])));
                    }

                    else if (i==test-1 && j==test-1 && k== test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(1,1,1)));
                    }

                    else if(i==test-1 && j==test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(1,1,tab_z[k])));
                    }

                    else if(k==test-1 && j==test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],1,1)));
                    }

                    else if(i==test-1 && k==test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(1,tab_y[j],1)));
                    }

                    else if(i==test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(1,tab_y[j],tab_z[k])));
                    }
                    else if(j==test-1 )
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],1,tab_z[k])));
                    }
                    else if(k==test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],tab_y[j],1)));
                    }
                    else
                        std::cout<<"impossible\n";

                    tmp_ligne_maj[j]=memoire;
                    tmp_box_maj[i][j]=etage_inf+memoire;

                    //comparer max
                    #pragma omp flush(discrepancytmp)
                    if(discrepancytmp<=tmp )
                    {
                    #pragma omp crititcal
                        if(discrepancytmp<=tmp  )
                        {
                            discrepancytmp=tmp;
                        }
                    }
				}///fin j
            }///fin i
        }///fin k

        /////////// BOUCLES SUPPLEMENTAIRES POUR TESTER LES BOITES VIDES SUR LES AXES X, Y ET Z///////////////////////////////

            #pragma omp parallel for num_threads(nb_threads)
            for(int i=0;i<pts_number;i++)
            {
                double tmp=box_area_calcul(tab_x[i],tab_y[0],tab_z[0]);
                #pragma omp flush(discrepancytmp)
                 if(discrepancytmp<=tmp)
                 {
                    #pragma omp crititcal
                    if(discrepancytmp<=tmp)
                    {
                        discrepancytmp=tmp;
                    }
                }
            }

            #pragma omp parallel for num_threads(nb_threads)
            for(int i=1;i<pts_number;i++)
            {
                double tmp=box_area_calcul(tab_x[0],tab_y[i],tab_z[0]);
                #pragma omp flush(discrepancytmp)
                 if(discrepancytmp<=tmp)
                 {
                    #pragma omp crititcal
                    if(discrepancytmp<=tmp)
                    {
                        discrepancytmp=tmp;
                    }
                }
            }

            #pragma omp parallel for num_threads(nb_threads)
            for(int i=1;i<pts_number;i++)
            {
                double tmp=box_area_calcul(tab_x[0],tab_y[0],tab_z[i]);
                #pragma omp flush(discrepancytmp)
                 if(discrepancytmp<=tmp)
                 {
                    #pragma omp crititcal
                    if(discrepancytmp<=tmp)
                    {
                        discrepancytmp=tmp;
                    }
                }
            }

        delete[] tmp_x;
        delete[] tmp_y;
        delete[] tmp_z;

        ////////////////////// FIN DU CALCUL //////////////////////////////

        return discrepancytmp;
}

void export_data_3D(double discrepancy){

            std::ofstream file(fn_output.c_str());
            if(!file)
            std::cerr << "Can't open ouptut" << std::endl;

            file << pts_number << "  " << discrepancy << std::endl;
            file.close();
}

/////////////////////////  SEQUENCE  ///////////////////////////////

bool sequence_import_data_2D(int palier){

    std::string line;

	std::ifstream fichier(fn_input.c_str(), std::ios::in);
    if(!fichier.fail())
	{
        std::cout << "  ***************************************  \n\n" ;

		std::cout << " Chargement des donnees...\n" ;
		pts_number = 0;

		while (std::getline(fichier, line)){
            if(line != "" && line !="\n" && line!=" "){
                pts_number++;
        }
        else continue;
	}

		std::cout <<" "<<pts_number<<" lignes. \n";
		if(pts_number>0){

			fichier.clear();
			fichier.seekg(0, std::ios::beg);

			double x=0;
			double y=0;
			int i=0;

            points_tab_2D = new Points_2D[pts_number];
            tab_x=new double[palier];
			tab_y=new double[palier];
			pos_y_point_axe_x=new double[palier];
			pos_x_point_axe_y=new double[palier];

			tmp_ligne=new int[palier];
			tmp_ligne_maj=new int[palier];

            while(fichier >> x >> y)
            {
                points_tab_2D[i].set_pos_x(x);
				points_tab_2D[i].set_pos_y(y);

				i++;
            }
			std::cout << "\n Data Ok !\n\n" ;

        }

			fichier.close();
			return true;
    }
    else {
        std::cerr << " Can't open input !\n" ;
        return false;
    }

}

void sequence_export_data(int i,double discrepancy){

        std::ofstream file(fn_output.c_str(),std::ios::app);
        if(!file)
        std::cerr << "Can't open ouptut" << std::endl;

        file << pts_number << "  " << discrepancy << std::endl;
        file.close();
}

double sequence_calcul_discrepancy_2D(int palier){

            double discrepancy=0;
////////////////////////// INITIALISATION TAB POS//////////////////////////

            for(int i=0;i<palier;i++)
            {
                #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<palier;j++)
                {
                    if(std::fabs(tab_x[i]-points_tab_2D[j].get_pos_x())<std::numeric_limits<double>::epsilon())
                        pos_y_point_axe_x[i]=points_tab_2D[j].get_pos_y(); //on a la position a partir duquel on trouve le point

                    if(std::fabs(tab_y[i]-points_tab_2D[j].get_pos_y())<std::numeric_limits<double>::epsilon())
                        pos_x_point_axe_y[i]=points_tab_2D[j].get_pos_x();
                }
            }

///////////////////////////////////////////////////////////////////////////

            #pragma omp parallel for num_threads(nb_threads)
			for(int i=0;i<palier;i++)
			{
                tmp_ligne[i]=0;
                tmp_ligne_maj[i]=0;
			}

            const int test = palier;

////////////////////////// DISCREPANCE //////////////////////////////////////

			for(int i=0;i<palier;i++)
			{

            #pragma omp parallel for num_threads(nb_threads)
            for(int k=0;k<palier;k++)
            {
                tmp_ligne[k]=tmp_ligne_maj[k];
            }

            #pragma omp parallel for firstprivate(tmp_ligne,pos_y_point_axe_x,pos_x_point_axe_y,tab_x,tab_y) num_threads(nb_threads)
            for(int j=0;j<test;j++)
            {
                    int memoire=0;

                    if(j>0)
                    memoire=tmp_ligne[j-1];

                    if(tab_y[j]>pos_y_point_axe_x[i])
                        {
                            memoire++;
                        }

                    if(tab_x[i]>pos_x_point_axe_y[j])
                        {
                            memoire++;
                        }

                    if(std::fabs(tab_x[i]-pos_x_point_axe_y[j])<std::numeric_limits<double>::epsilon() && std::fabs(tab_y[j]-pos_y_point_axe_x[i])<std::numeric_limits<double>::epsilon())
                        {
                            memoire++;
                        }

                    double tmp=0;

                    if(i<test-1 && j<test-1)
                    {
                        tmp=fabs((((double)memoire))/((double)test)-box_area_calcul(tab_x[i],tab_y[j]));
                        tmp=std::max(tmp,fabs((((double)memoire))/((double)test)-box_area_calcul(tab_x[i+1],tab_y[j+1])));
                        tmp=std::max(tmp,fabs((((double)memoire))/((double)test)-box_area_calcul(tab_x[i+1],tab_y[j])));
                        tmp=std::max(tmp,fabs((((double)memoire))/((double)test)-box_area_calcul(tab_x[i],tab_y[j+1])));
                    }

                    else if (i==test-1 && j==test-1)
                        tmp=std::max(tmp,fabs((((double)memoire))/((double)test)-box_area_calcul(1,1)));

                    else if(i==test-1)
                        tmp=std::max(tmp,fabs((((double)memoire))/((double)test)-box_area_calcul(1,tab_y[j])));

                    else if (j==test-1)
                        tmp=std::max(tmp,fabs((((double)memoire))/((double)test)-box_area_calcul(tab_x[i],1)));

                    else
                        std::cout<<"impossible\n";

                    tmp_ligne_maj[j]=memoire;

                    #pragma omp flush(discrepancy)
                    if(discrepancy<=tmp)
                    {
                    #pragma omp crititcal
                        if(discrepancy<=tmp)
                        {
                            discrepancy=tmp;
                        }
                    }
				}
            }


/////////// BOUCLES SUPPLEMENTAIRES POUR TESTER LES BOITES VIDES SUR LES AXES X ET Y///////////////////////////////

            #pragma omp parallel for num_threads(nb_threads)
            for(int i=0;i<palier;i++)
            {
                double tmp=box_area_calcul(tab_x[i],tab_y[0]);
                #pragma omp flush(discrepancy)
                 if(discrepancy<=tmp)
                 {
                    #pragma omp crititcal
                    if(discrepancy<=tmp)
                    {
                        discrepancy=tmp;
                    }
                }
            }

            #pragma omp parallel for num_threads(nb_threads)
            for(int i=1;i<palier;i++)
            {
                double tmp=box_area_calcul(tab_x[0],tab_y[i]);
                #pragma omp flush(discrepancy)
                 if(discrepancy<=tmp)
                 {
                    #pragma omp crititcal
                    if(discrepancy<=tmp)
                    {
                        discrepancy=tmp;
                    }
                }
            }

////////////////////// FIN DU CALCUL //////////////////////////////

    return discrepancy;
}

bool sequence_import_data_1D(int palier){

    std::string line;
	std::ifstream fichier(fn_input.c_str(), std::ios::in);
    if(!fichier.fail())
	{
        std::cout << "  ***************************************  \n\n" ;
		std::cout << " Chargement des donnees...\n" ;
		pts_number = 0;

		while (std::getline(fichier, line)){
            if(line != "" && line !="\n" && line!=" "){
                pts_number++;
            }
            else continue;
        }

		std::cout <<" "<<pts_number<<" lignes. \n";
		if(pts_number>0){

			fichier.clear();
			fichier.seekg(0, std::ios::beg);

			double x=0;
			int i=0;

			points_tab_1D = new Points_1D[pts_number];
			tab_x=new double[palier];

            while(fichier >> x )
            {
                points_tab_1D[i].set_pos_x(x);
				i++;
            }
			std::cout << "\n Data Ok !\n\n" ;
        }

        fichier.close();
        return true;
    }
    else {
        std::cerr << " Can't open input !\n" ;
        return false;
    }

}

double sequence_calcul_discrepancy_1D(int palier){

            double discrepancy=0;
			quickSort(tab_x,0,palier-1);
            const int test = palier;

////////////////////////// DISCREPANCE //////////////////////////////

            #pragma omp parallel for num_threads(nb_threads)
            for(int j=0;j<test;j++)
            {
                int memoire=j;
                double tmp=fabs(tab_x[j]-((2*(j+1)-1)/(double)(2*palier)));

                #pragma omp flush(discrepancy)
                if(discrepancy<=tmp)
                {
                    #pragma omp crititcal
                    if(discrepancy<=tmp)
                    {
                        discrepancy=tmp;
                    }
                }
            }

            discrepancy=discrepancy+((double)1/(double)(2*palier));

////////////////////// FIN DU CALCUL //////////////////////////////

           return discrepancy;
}

bool sequence_import_data_3D(int palier){

        std::string line;

        std::ifstream fichier(fn_input.c_str(), std::ios::in);
    	if(!fichier.fail())
        {
            std::cout << "  ***************************************  \n\n" ;
            std::cout << " Chargement des donnees...\n" ;
            pts_number = 0;

            while (std::getline(fichier, line)){
                if(line != "" && line !="\n" && line!=" "){
                    pts_number++;
                }
            else continue;
            }

            std::cout <<" "<<pts_number<<" lignes. \n";
            if(pts_number>0){

                fichier.clear();
                fichier.seekg(0, std::ios::beg);

                double x=0;
                double y=0;
                double z=0;
                int i=0;

                points_tab_3D = new Points_3D[pts_number];
                tab_x=new double[palier];
                tab_y=new double[palier];
                tab_z=new double[palier];

                pos_y_point_axe_x=new double[palier];
                pos_y_point_axe_z=new double[palier];

                pos_x_point_axe_y=new double[palier];
                pos_x_point_axe_z=new double[palier];

                pos_z_point_axe_x=new double[palier];
                pos_z_point_axe_y=new double[palier];

                tmp_ligne=new int[palier];
                tmp_ligne_maj=new int[palier];

                tmp_box=new int*[palier];
                tmp_box_maj=new int*[palier];

                for(int i=0;i<palier;i++)
                {
                    tmp_box[i]=new int [palier];
                    tmp_box_maj[i]=new int [palier];

                }

                while(fichier >> x >> y >> z)
                {
                    points_tab_3D[ i ].set_pos_x( x );
                    points_tab_3D[ i ].set_pos_y( y );
                    points_tab_3D[ i ].set_pos_z( z );


                    i++;
                }

                std::cout << "\n Data Ok !\n\n" ;

            }

			fichier.close();
			return true;
    }
    else {
        std::cerr << " Can't open input !\n" ;
        return false;
    }

}

double sequence_calcul_discrepancy_3D(int palier){

            double discrepancytmp=0;

////////////////////////// INITIALISATION TAB POS///////////////////////////////

            bool x,y,z;
            double* tmp_x = new double[palier];
            double* tmp_y = new double[palier];
            double* tmp_z = new double[palier];

            for(int i=0;i<palier;i++)
            {
                tmp_x[i]=0;
                tmp_y[i]=0;
                tmp_z[i]=0;
            }


           for(int i=0;i<palier;i++)
            {
                x=true;
                y=true;
                z=true;
                #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<palier;j++)
                {
                    if(x==true){
                        if(std::fabs(tab_x[i]-points_tab_3D[j].get_pos_x())<std::numeric_limits<double>::epsilon() && tmp_x[j]==0)
                        {
                            pos_y_point_axe_x[i]=points_tab_3D[j].get_pos_y();
                            pos_z_point_axe_x[i]=points_tab_3D[j].get_pos_z();
                            x=false;
                            tmp_x[j]=5;
                        }
                    }
                    if(y==true){
                        if(std::fabs(tab_y[i]-points_tab_3D[j].get_pos_y())<std::numeric_limits<double>::epsilon() && tmp_y[j]==0)
                        {
                            pos_x_point_axe_y[i]=points_tab_3D[j].get_pos_x();
                            pos_z_point_axe_y[i]=points_tab_3D[j].get_pos_z();
                            y=false;
                            tmp_y[j]=5;
                        }
                    }
                    if(z==true){
                        if(std::fabs(tab_z[i]-points_tab_3D[j].get_pos_z())<std::numeric_limits<double>::epsilon() && tmp_z[j]==0)
                        {
                            pos_x_point_axe_z[i]=points_tab_3D[j].get_pos_x();
                            pos_y_point_axe_z[i]=points_tab_3D[j].get_pos_y();
                            z=false;
                            tmp_z[j]=5;
                        }
                    }
                }
            }

///////////////////////////////////////////////////////////////////////////

            #pragma omp parallel for num_threads(nb_threads)
			for(int i=0;i<palier;i++)
			{
                tmp_ligne[i]=0;
                tmp_ligne_maj[i]=0;
			}
/////////////////////////////////////////////////////////////////////////////////////////

            for(int i=0;i<palier;i++)
            {
              #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<palier;j++)
                {
                    tmp_box[i][j]=0;
                    tmp_box_maj[i][j]=0;
                }
            }

            const int test = palier;
////////////////////////// DISCREPANCE //////////////////////////////////////////////////

    for(int k=0; k<palier;k++){

            for(int o=0; o<palier;o++){
                #pragma omp parallel for
                for(int p=0; p<palier;p++){
                    tmp_box[o][p]=tmp_box_maj[o][p];

                }
            }
                #pragma omp parallel for num_threads(nb_threads)
                for(int l=0;l<palier;l++)
                {
                    tmp_ligne_maj[l]=0;
                }

			for(int i=0;i<palier;i++)
			{
                #pragma omp parallel for num_threads(nb_threads)
                for(int l=0;l<palier;l++)
                {
                    tmp_ligne[l]=tmp_ligne_maj[l];
                }


                #pragma omp parallel for firstprivate(tmp_ligne,tmp_box,pos_z_point_axe_x,pos_z_point_axe_y,pos_y_point_axe_x,pos_y_point_axe_z,pos_x_point_axe_y,pos_x_point_axe_z,tab_x,tab_y,tab_z) num_threads(nb_threads)
				for(int j=0;j<test;j++)
				{

                   int memoire=0;
                   int etage_inf=0;

                    if(j>0 )
                    {
                        memoire=tmp_ligne[j-1];
                    }

                    etage_inf=tmp_box[i][j];

                    if(tab_y[j]>pos_y_point_axe_x[i] && std::fabs(tab_z[k]-pos_z_point_axe_x[i])<std::numeric_limits<double>::epsilon())//////////////////////////////
                        {
                            memoire++;
                         }

                    if(tab_x[i]>pos_x_point_axe_y[j] && std::fabs(tab_z[k]-pos_z_point_axe_y[j])<std::numeric_limits<double>::epsilon())//////////////////////////////
                        {
                            memoire++;
                         }
                //si il y a un point dans la boîte
                if(
                    std::fabs(tab_x[i]-pos_x_point_axe_y[j])<std::numeric_limits<double>::epsilon() &&
                    std::fabs(tab_x[i]-pos_x_point_axe_z[k])<std::numeric_limits<double>::epsilon() &&

                    std::fabs(tab_y[j]-pos_y_point_axe_x[i])<std::numeric_limits<double>::epsilon() &&
                    std::fabs(tab_y[j]-pos_y_point_axe_z[k])<std::numeric_limits<double>::epsilon() &&

                    std::fabs(tab_z[k]-pos_z_point_axe_x[i])<std::numeric_limits<double>::epsilon() &&
                    std::fabs(tab_z[k]-pos_z_point_axe_y[j])<std::numeric_limits<double>::epsilon()
                )
                {
                    memoire++;
                }
                    double tmp=0;

                    tmp=fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],tab_y[j],tab_z[k]));
                //on check les autres angles de la boites
                   if(i<test-1 && j<test-1 && k<test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i+1],tab_y[j],tab_z[k])));
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],tab_y[j+1],tab_z[k])));
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i+1],tab_y[j+1],tab_z[k])));
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],tab_y[j],tab_z[k+1])));
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i+1],tab_y[j],tab_z[k+1])));
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],tab_y[j+1],tab_z[k+1])));
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i+1],tab_y[j+1],tab_z[k+1])));
                    }

                    else if (i==test-1 && j==test-1 && k== test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(1,1,1)));
                    }

                    else if(i==test-1 && j==test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(1,1,tab_z[k])));
                    }

                    else if(k==test-1 && j==test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],1,1)));
                    }

                    else if(i==test-1 && k==test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(1,tab_y[j],1)));
                    }

                    else if(i==test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(1,tab_y[j],tab_z[k])));
                    }
                    else if(j==test-1 )
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],1,tab_z[k])));
                    }
                    else if(k==test-1)
                    {
                        tmp=std::max(tmp,fabs((((double)etage_inf+memoire))/((double)test)-box_area_calcul(tab_x[i],tab_y[j],1)));
                    }
                    else
                        std::cout<<"impossible\n";

                    tmp_ligne_maj[j]=memoire;
                    tmp_box_maj[i][j]=etage_inf+memoire;

                    //comparer max
                    #pragma omp flush(discrepancytmp)
                    if(discrepancytmp<=tmp )
                    {
                    #pragma omp crititcal
                        if(discrepancytmp<=tmp  )
                        {
                            discrepancytmp=tmp;
                        }
                    }
				}///fin j
            }///fin i
        }///fin k

        /////////// BOUCLES SUPPLEMENTAIRES POUR TESTER LES BOITES VIDES SUR LES AXES X, Y ET Z///////////////////////////////

            #pragma omp parallel for num_threads(nb_threads)
            for(int i=0;i<palier;i++)
            {
                double tmp=box_area_calcul(tab_x[i],tab_y[0],tab_z[0]);
                #pragma omp flush(discrepancytmp)
                 if(discrepancytmp<=tmp)
                 {
                    #pragma omp crititcal
                    if(discrepancytmp<=tmp)
                    {
                        discrepancytmp=tmp;
                    }
                }
            }

            #pragma omp parallel for num_threads(nb_threads)
            for(int i=1;i<palier;i++)
            {
                double tmp=box_area_calcul(tab_x[0],tab_y[i],tab_z[0]);
                #pragma omp flush(discrepancytmp)
                 if(discrepancytmp<=tmp)
                 {
                    #pragma omp crititcal
                    if(discrepancytmp<=tmp)
                    {
                        discrepancytmp=tmp;
                    }
                }
            }

            #pragma omp parallel for num_threads(nb_threads)
            for(int i=1;i<palier;i++)
            {
                double tmp=box_area_calcul(tab_x[0],tab_y[0],tab_z[i]);
                #pragma omp flush(discrepancytmp)
                 if(discrepancytmp<=tmp)
                 {
                    #pragma omp crititcal
                    if(discrepancytmp<=tmp)
                    {
                        discrepancytmp=tmp;
                    }
                }
            }

        delete[] tmp_x;
        delete[] tmp_y;
        delete[] tmp_z;

////////////////////// FIN DU CALCUL //////////////////////////////

        return discrepancytmp;
}

/////////////////////////  MAIN  /////////////////////////////////

int main(int argc, char** argv){

double discrepancy=0;

import_options(argc,argv);

srand(time(NULL));

td=time(NULL);

if(!sequence){

    if(dimension == 1){

        if(import_data_1D()){

                discrepancy=calcul_discrepancy_1D();

                std::cout<<" Discrepancy : "<<discrepancy<<"\n";

                export_data_1D(discrepancy);

            }
    }

    else if(dimension == 2){

        if(import_data_2D()){


                discrepancy=calcul_discrepancy_2D();

                std::cout<<" Discrepancy : "<<discrepancy<<"\n";

                export_data_2D(discrepancy);

            }
    }

    else if (dimension == 3){

    if(import_data_3D()){

			discrepancy=calcul_discrepancy_3D();

            std::cout<<" Discrepancy : "<<discrepancy<<"\n";

            export_data_3D(discrepancy);

            }
    }

    else std::cout<<" Dimension not yet implemented. \n";

	 tf=time(NULL);

    std::cout << " Fin de la recherche en " << difftime(tf, td) << " secondes.\n" ;
}

else {

    if(dimension == 1){

        if(sequence_import_data_1D(1)){

            std::cout<<"Sequence : do stuff \n";

            std::ofstream fichier(fn_output.c_str(),std::ios::out);
            if(fichier)
            {
                // fichier<<"Discrepance pour "<<pts_number<<" points.\n";
                fichier.close();
            }
            else
                std::cerr<<"Erreur à l'ouverture fichier output !"<<std::endl;

            Points_1D* sequence_points_tab_1D = new Points_1D[pts_number];
            int mem_x = -1;
            #pragma omp parallel for num_threads(nb_threads)
            for(int i=0;i<pts_number;i++)
            {
                sequence_points_tab_1D[i]=points_tab_1D[i];
            }

            delete [] points_tab_1D;

            for(int palier=1;palier<=pts_number;palier++)
            {
                double* sequence_tab_x;

                if(palier>1){
                    sequence_tab_x=new double[palier-1];
                        #pragma omp parallel for num_threads(nb_threads)
                        for(int j=0;j<palier-1;j++)
                        {
                            sequence_tab_x[j]=tab_x[j];
                        }

                }

            bool x_insert;
            x_insert=false;
            tab_x=new double[palier];
            points_tab_1D=new Points_1D[palier];



            if(palier==1)
            {
                points_tab_1D[0].set_pos_x(sequence_points_tab_1D[0].get_pos_x());
                tab_x[0]=sequence_points_tab_1D[0].get_pos_x();
            }

            else{

                #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<palier;j++)
                {
                    points_tab_1D[j].set_pos_x(sequence_points_tab_1D[j].get_pos_x());
                }
                #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<palier-1;j++)
                {
                    if(sequence_tab_x[j]<sequence_points_tab_1D[palier-1].get_pos_x())
                    {
                        tab_x[j]=sequence_tab_x[j];
                    }
                    else if(!x_insert && sequence_tab_x[j]>sequence_points_tab_1D[palier-1].get_pos_x() )
                    {
                        x_insert=true;
                        mem_x=j;
                    }
                }



                if(mem_x!=-1)
                {
                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=0;j<mem_x;j++)
                    {
                        tab_x[j]=sequence_tab_x[j];
                    }

                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=mem_x+1;j<palier;j++)
                    {
                        tab_x[j]=sequence_tab_x[j-1];
                    }
                    tab_x[mem_x]=sequence_points_tab_1D[palier-1].get_pos_x();
                }
                else
                    tab_x[palier-1]=sequence_points_tab_1D[palier-1].get_pos_x();

            }
            discrepancy= sequence_calcul_discrepancy_1D(palier);
            sequence_export_data(palier,discrepancy);
            }
        }
    }

    else if(dimension == 2){

        if(sequence_import_data_2D(1)){

        std::cout<<"Sequence : do stuff \n";

            std::ofstream fichier(fn_output.c_str(),std::ios::out);
            if(fichier)
            {
                // fichier<<"Discrepance pour "<<pts_number<<" points.\n";
                fichier.close();
            }
            else
                std::cerr<<"Erreur à l'ouverture fichier output !"<<std::endl;

        Points_2D* sequence_points_tab_2D = new Points_2D[pts_number];

        #pragma omp parallel for num_threads(nb_threads)
        for(int i=0;i<pts_number;i++)
        {
            sequence_points_tab_2D[i]=points_tab_2D[i];
        }

        delete [] points_tab_2D;

        for(int palier=1;palier<=pts_number;palier++)
        {
/////////////////////////////////////////////////////////////
            double* sequence_tab_x;
            double* sequence_tab_y;

            if(palier>1){
                sequence_tab_x=new double[palier-1];
                sequence_tab_y=new double[palier-1];

                #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<palier-1;j++)
                {
                    sequence_tab_x[j]=tab_x[j];
                    sequence_tab_y[j]=tab_y[j];
                }
            }

            tab_x=new double[palier];
			tab_y=new double[palier];

			pos_y_point_axe_x=new double[palier];
			pos_x_point_axe_y=new double[palier];
			tmp_ligne=new int[palier];
			tmp_ligne_maj=new int[palier];

            points_tab_2D=new Points_2D[palier];

            bool x_insert,y_insert;
            x_insert=false;
            y_insert=false;
            int mem_x = -1;
            int mem_y = -1;

//////////////////////////////////////////////////////////////


            if(palier==1)
            {
                points_tab_2D[0].set_pos_x(sequence_points_tab_2D[0].get_pos_x());
                points_tab_2D[0].set_pos_y(sequence_points_tab_2D[0].get_pos_y());

                tab_x[0]=sequence_points_tab_2D[0].get_pos_x();
                tab_y[0]=sequence_points_tab_2D[0].get_pos_y();
            }

            else{
                #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<palier;j++)//on ajoute le nv point ds le tab de points
                {
                    points_tab_2D[j].set_pos_x(sequence_points_tab_2D[j].get_pos_x());
                    points_tab_2D[j].set_pos_y(sequence_points_tab_2D[j].get_pos_y());
                }

                #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<palier-1;j++)//on remplit ceux d'avant
                {
                    if(sequence_tab_x[j]<sequence_points_tab_2D[palier-1].get_pos_x())
                    {
                        tab_x[j]=sequence_tab_x[j];
                    }
                    else if(!x_insert && sequence_tab_x[j]>sequence_points_tab_2D[palier-1].get_pos_x() )
                    {
                        x_insert=true;
                        mem_x=j;
                    }


                    if(sequence_tab_y[j]<sequence_points_tab_2D[palier-1].get_pos_y())
                    {
                        tab_y[j]=sequence_tab_y[j];
                    }
                    else if(!y_insert && sequence_tab_y[j]>sequence_points_tab_2D[palier-1].get_pos_y() )
                    {
                        y_insert=true;
                        mem_y=j;
                    }
                }


                if(mem_x!=-1)
                {
                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=0;j<mem_x;j++)
                    {
                        tab_x[j]=sequence_tab_x[j];
                    }

                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=mem_x+1;j<palier;j++)
                    {
                        tab_x[j]=sequence_tab_x[j-1];
                    }


                    tab_x[mem_x]=sequence_points_tab_2D[palier-1].get_pos_x();
                }
                else
                    tab_x[palier-1]=sequence_points_tab_2D[palier-1].get_pos_x();



                if(mem_y!=-1)
                {
                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=0;j<mem_y;j++)
                    {
                        tab_y[j]=sequence_tab_y[j];
                    }

                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=mem_y+1;j<palier;j++)
                    {
                        tab_y[j]=sequence_tab_y[j-1];
                    }


                    tab_y[mem_y]=sequence_points_tab_2D[palier-1].get_pos_y();
                    }
                    else
                        tab_y[palier-1]=sequence_points_tab_2D[palier-1].get_pos_y();


            }
            discrepancy= sequence_calcul_discrepancy_2D(palier);
            sequence_export_data(palier,discrepancy);
        }

    }
}

    else if (dimension == 3){

    if(sequence_import_data_3D(1)){

        std::cout<<"Sequence : do stuff \n";

        std::ofstream fichier(fn_output.c_str(),std::ios::out);
        if(fichier)
        {
            // fichier<<"Discrepance pour "<<pts_number<<" points.\n";
            fichier.close();
        }
        else
            std::cerr<<"Erreur à l'ouverture fichier output !"<<std::endl;

        Points_3D* sequence_points_tab_3D = new Points_3D[pts_number];

        #pragma omp parallel for num_threads(nb_threads)
        for(int i=0;i<pts_number;i++)
        {
            sequence_points_tab_3D[i]=points_tab_3D[i];
        }

        delete [] points_tab_3D;

        for(int palier=1;palier<=pts_number;palier++)
        {
            double* sequence_tab_x;
            double* sequence_tab_y;
            double* sequence_tab_z;

            if(palier>1){
                sequence_tab_x=new double[palier-1];
                sequence_tab_y=new double[palier-1];
                sequence_tab_z=new double[palier-1];

                #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<palier-1;j++)
                {
                    sequence_tab_x[j]=tab_x[j];
                    sequence_tab_y[j]=tab_y[j];
                    sequence_tab_z[j]=tab_z[j];
                }
            }

            tab_x=new double[palier];
			tab_y=new double[palier];
			tab_z=new double[palier];

			pos_x_point_axe_y=new double[palier];
			pos_x_point_axe_z=new double[palier];

			pos_y_point_axe_x=new double[palier];
			pos_y_point_axe_z=new double[palier];

			pos_z_point_axe_x=new double[palier];
			pos_z_point_axe_y=new double[palier];

			tmp_ligne=new int[palier];
			tmp_ligne_maj=new int[palier];

            tmp_box=new int*[palier];
            tmp_box_maj=new int*[palier];

            #pragma omp parallel for num_threads(nb_threads)
            for(int i=0;i<palier;i++)
            {
                    tmp_box[i]=new int [palier];
                    tmp_box_maj[i]=new int [palier];

            }

             points_tab_3D=new Points_3D[palier];

             bool x_insert,y_insert,z_insert;
             x_insert=false;
             y_insert=false;
             z_insert=false;
             int mem_x = -1;
             int mem_y = -1;
             int mem_z = -1;

//////////////////////////////////////////////////////////////

            if(palier==1)
            {
                points_tab_3D[0].set_pos_x(sequence_points_tab_3D[0].get_pos_x());
                points_tab_3D[0].set_pos_y(sequence_points_tab_3D[0].get_pos_y());
                points_tab_3D[0].set_pos_z(sequence_points_tab_3D[0].get_pos_z());

                tab_x[0]=sequence_points_tab_3D[0].get_pos_x();
                tab_y[0]=sequence_points_tab_3D[0].get_pos_y();
                tab_z[0]=sequence_points_tab_3D[0].get_pos_z();
            }

            else{

                #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<palier;j++)
                {
                    points_tab_3D[j].set_pos_x(sequence_points_tab_3D[j].get_pos_x());
                    points_tab_3D[j].set_pos_y(sequence_points_tab_3D[j].get_pos_y());
                    points_tab_3D[j].set_pos_z(sequence_points_tab_3D[j].get_pos_z());
                }

                #pragma omp parallel for num_threads(nb_threads)
                for(int j=0;j<palier-1;j++)
                {
                    if(sequence_tab_x[j]<sequence_points_tab_3D[palier-1].get_pos_x())
                    {
                        tab_x[j]=sequence_tab_x[j];
                    }
                    else if(!x_insert && sequence_tab_x[j]>sequence_points_tab_3D[palier-1].get_pos_x() )
                    {
                        x_insert=true;
                        mem_x=j;
                    }


                    if(sequence_tab_y[j]<sequence_points_tab_3D[palier-1].get_pos_y())
                    {
                        tab_y[j]=sequence_tab_y[j];
                    }
                    else if(!y_insert && sequence_tab_y[j]>sequence_points_tab_3D[palier-1].get_pos_y() )
                    {
                        y_insert=true;
                        mem_y=j;
                    }


                    if(sequence_tab_z[j]<sequence_points_tab_3D[palier-1].get_pos_z())
                    {
                        tab_z[j]=sequence_tab_z[j];
                    }
                    else if(!z_insert && sequence_tab_z[j]>sequence_points_tab_3D[palier-1].get_pos_z() )
                    {
                        z_insert=true;
                        mem_z=j;
                    }
                }



                if(mem_x!=-1)
                {
                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=0;j<mem_x;j++)
                    {
                        tab_x[j]=sequence_tab_x[j];
                    }

                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=mem_x+1;j<palier;j++)
                    {
                        tab_x[j]=sequence_tab_x[j-1];
                    }
                    tab_x[mem_x]=sequence_points_tab_3D[palier-1].get_pos_x();
                }
                else
                    tab_x[palier-1]=sequence_points_tab_3D[palier-1].get_pos_x();



                if(mem_y!=-1)
                {
                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=0;j<mem_y;j++)
                    {
                        tab_y[j]=sequence_tab_y[j];
                    }

                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=mem_y+1;j<palier;j++)
                    {
                        tab_y[j]=sequence_tab_y[j-1];
                    }
                    tab_y[mem_y]=sequence_points_tab_3D[palier-1].get_pos_y();
                    }
                    else
                        tab_y[palier-1]=sequence_points_tab_3D[palier-1].get_pos_y();



                if(mem_z!=-1)
                {
                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=0;j<mem_z;j++)
                    {
                        tab_z[j]=sequence_tab_z[j];
                    }

                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=mem_z+1;j<palier;j++)
                    {
                        tab_z[j]=sequence_tab_z[j-1];
                    }
                    tab_z[mem_z]=sequence_points_tab_3D[palier-1].get_pos_z();
                    }
                    else
                        tab_z[palier-1]=sequence_points_tab_3D[palier-1].get_pos_z();
                }
                discrepancy= sequence_calcul_discrepancy_3D(palier);
                sequence_export_data(palier,discrepancy);

              for(int i=0;i<palier;i++)
                {
                    delete[] tmp_box[i];
                    delete[] tmp_box_maj[i];

                }
            }
        }
    }


}

}











