# PESA

Paramètres spécifiques:
	* Taille de l'archive
	* Taille d'une case de l'hyper-grille

À chaque itération:
	* Récupérer les individus non dominés de la population
	* Sinon, pour chacun de ces individus
		* Si il est dominé par un individu de l'archive, passer au suivant
		* Sinon
			* L'ajouter à l'archive
			* Retirer les individus de l'archive dominés par ce nouvel individu
			* Si l'archive est pleine, on retire l'un des individus de l'archive avec le plus grand squeeze factor
	* Si critère d'arrêt, retourner l'archive
	* On génère une nouvelle population
	* Tant que la population n'est pas complète
		* Effectuer un crossover avec une probabilité Pc (en sélectionnant deux individus en fonction de leur squeeze factor)
		* Sinon, effectuer une mutation (en sélectionnant un individu de la population)

   
# IBEA
Paramètres spécifiques:
	* kappa
	(* l'indicateur binaire, fixé pour nous)

À chaque itération:
	* Calculer la fitness de chaque individu en fonction de l'indicateur binaire
	* Si la population compte plus que N éléments
		* Récupérer le pire individu et le retirer (ne pas oublier de mettre à jour la fitness des individus restants)
	* Si critère d'arrêt, retourner les individus non dominés
	* On génère une nouvelle population (crossover et mutation)
