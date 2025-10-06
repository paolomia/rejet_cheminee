# Rapport interne - Calcul de conduits de fumée pour AF005215

# 1. Introduction


## 1.1 Hypothèses
Le texte de réference est la norme NF EN 13384-1+A1:2019 "Conduits de fumée - Métho
des de calcul"


Le système prise en examen comporte :

- Un bruleur gaz à air soufflé Riello BS3/M
- Un foyer échangeur de fabrication Omia GQ01 63
- Un conduit de raccordement de 200 mm de long et 200 mm de diamètre avec un angle de 90°
- un conduit de fumée de longuer 6.3 m et de diamètre 200 mm à l'intérieur du bâtiment
- un conduit de fumée de longuer 1.2 m et de diamètre 200 mm à l'extérieur du bâtiment
- L'hauteur totale du conduit de fumée est donc de 7.5 m

Les conduits sont en tôle métallique non isolé, 1 mm d'épaisseur. Le conduit situé à l'intérieur du batimen est dans un environnement chauffé (15 °C).

Aucun volet n'est installé dans le conduit de fumée.
La puissance du bruleur n'est pas modulable, le bruleur fonctionne en tout ou rien.
Un cône de condensation est situé au fond du raccord, avec un trou pour l'évacuation des condensats.

L'adresse de l'installation est Havighorster Weg 8C, 21031 Hamburg, Germania. Le site est situé à plus de 20 km de la mer et à une altitude de 41 m au-dessus du niveau de la mer.
Les températures minime des dernières années sont les suivantes :
2025 	-9.9 °C
2024 	-9.1 °C
2023 	-8.5 °C
2022 	-9.1 °C
2021 	-14.4 °C

La sortie de fumée n'est pas dans une zone soumise a une pression défavorable (5.10.4 pg. 36 de la norme).

D'après les plans, la puissance nominale de l'installation est de 140 kW
![img_1.png](img_1.png)

D'après essais réalisés par Omia, le rendement est de 86 % et la température des fumées à la sortie du foyer est de 310 °C.
![img_2.png](img_2.png)

Je n'ai pas de données concernant la pression maximale supportée par le conduit de fumée.

Nous n'avons pas de données sur le debit massique des fumées, mais la norme (5.5.2.1 pg. 20 et Tab B.1 pg. 99) donne une formule pour l'estimation théorique. Pour les valeurs données, nous avons un debit massique estimé de 0.068 kg/s.

Le chauffage marche en tout ou rien, donc nous pouvons prendre les valeurs de débit massique et de température des fumées qui correspondent à la puissance nominale.

## 1.2 Critères

Nous sommes dans le cas d'un conduit fonctionnant sous pression positive. Nous devons respecter les 4 critères suivants (pg.17, 5.1)

### 1.2 Critères de pression

1. La pression positive maximale à l’entrée des fumées dans le conduit doit être ≤ à la pression différentielle maximale à l’entrée des fumées dans le conduit.

Donc P_ZO ≤ P_ZOe 

avec :

Pression positive maximale (P_ZO) : pression que les fumées exercent à l’entrée du conduit de fumée à cause de l'effet du tirage, du vent et de pertes de charge du circuit de fumée.

Pression différentielle maximale (P_ZOe) : limite maximale de différence de pression que l’appareil de combustion peut accepter pour fonctionner correctement.

2. La pression positive maximale dans le conduit de raccordement des fumées et dans le conduit de fumée ne doit pas être supérieure à la pression maximale pour laquelle ils ont été désignés

Donc P_ZO ≤ P_Zexcess

avec :

P_Zexcess : pression maximale que le conduit de fumée peut supporter en service.

3. La pression positive minimale au niveau de l’admission des fumées dans le conduit doit être supérieure ou égale à la pression différentielle minimale au niveau de l’admission des fumées dans le conduit

Donc P_ZOmin ≥ P_ZOemin 

avec :

Pression positive minimale (PZOmin) : plus petite surpression que l’on peut avoir à l’entrée du conduit dans les conditions les plus défavorables (ex. extérieur très froid, pertes de charge plus fortes)

Pression différentielle minimale (P_ZOemin) : différence de pression minimale que l’appareil de combustion a besoin pour fonctionner correctement

### 1.2 Critères de température

4. La température de la paroi intérieure à la sortie du conduit de fumée doit être supérieure ou égale à la limite de température.
Donc T_iob ≥ T_ig

avec :

Température de paroi intérieure à la sortie (T_iob) : température réelle de la paroi interne du conduit de fumée à la sortie

Limite de température (T_ig) : c’est la température minimale fixée par la norme

Le valeur de T_ig est :
- égale à la température de rosée dans le cas d'un calcul sans condensation
- égale à 0 °C dans le cas d'un calcul avec condensation

## 1.3 Prises de mesure 03/10/2025

Avec Francis Diot, nous avons réalisé des prises de mesure le 03/10/2025 sur le foyer échangeur gaz de l'usine.

Température ambiante (non enregistrée) : ~25 °C

Pression dans la chambre de combustion pendant le fonctionnement : entre 0 et -10 Pa

Pression prise dans le conduit de fumée au niveau du trou de condensation: entre -30 et -40 Pa

Température de fumée au niveau du trou de condensation : 310 °C

C02 au niveau du trou de condensation : 8.5 %

Note:

- la pression dans la chambre de combustion monte considerablement lors de l'allumage du bruleur.
- D'après le test avec fumigène, la pression dans le conduit est négative même au bruleur éteint. La valeur n'est pas mésurable avec l'appareil de Francis.
- La difference de pression nous donne une estimation de la perte de charge dans le foyer échangeur, qui est d'environ 30 Pa.

## 1.4 Remarques

La norme prevoit des critères different pour le cas d'un système de chauffage à pression positive ou négative. J'hésite en partie à situer le notre, car je vois des arguments pour les deux cas.

**Pression positive:**
- Le bruleur est un bruleur à air soufflé, avec une plage de fonctionnnement en pression positive plus importante.
- Exigences moins strictes pour les calculs

**Pression négative:**
- Le bruleur fonctionne en pression négative
- Moins d'exigeances pour la ténue du conduit de fumée

J'ai choisi de faire les calculs en pression positive.

# 2. Calculs

## 2.1 Pression - Critères 1-3

Nous avons les données suivantes (voir note de calcul):

P_ZO ≈ -31 Pa
P_ZOmin ≈ -43 Pa


La perte de charge du foyer échangeur est estimée à 30 Pa. Étant donné l'incertitude sur cette valeur, nous pouvons prendre une marge de securité du 100% et affirmer que la perte de charge du foyer échangeur est comprise entre 15 et 60 Pa.

Pour les plages de pression du bruleur Riello (modèle BS3/M), nous avons les données suivantes (extrait de la documentation technique):
![img.png](img.png)

donc entre -50 Pa et 200 Pa à la puissance nominale.

Donc, la plage de valeurs de pressions admissibles au niveau de l'entrée du conduit de fumée est :

- Valeur minimale : La plus petite valeur de pression admise pour le bruleurs moins la plus petite valeur de perte de charge du foyer échangeur.
P_ZOemin = -50 - 15 = -65 Pa

- Valeur maximale : La plus grande valeur de pression pour le bruleurs moins la plus grande valeur de perte de charge du foyer échangeur.
P_ZOe = 200 - 60 = 140 Pa

Nous avons donc :
- P_ZO ≈ -31 Pa ≤ P_ZOe = 140 Pa  (**Critère 1 : OK**)
- P_ZO_min = -43 Pa ≥ P_ZOemin = -65 Pa (**Critère 3 : OK**)


## 2.2 Pression - Critère 2

Je n'ai pas d'élément concernant la pression maximale supportée par nos conduits de fumée. 

La presence du trou pour l'évacuation des condensats dans le raccordement implique que la pression doit être négative pour assurer qu'il n'y ait pas de sortie de fumée. Étant donné que la pression au niveau du raccordement est comprise entre -26 et -38 Pa (la perte de charge du raccordement est de 5 Pa) les conditions de pression sont théoriquement respectées.

Dans le cadre de la norme, l'action du vent (25 Pa) doit être prise en compte uniquement dans le cas d'une position de cheminée défavorable.
Sauf que le batiment en question est en depression, à cause de l'extraction des aires de préparation.

Mais cela nous engage à garantir que la pression dans le conduit de fumée reste toujours négative.

**Critère 2: ???**

## 2.3 Température - Critère 4

La température de rosée est estimée à 53 °C, la température de paroi intérieure à la sortie du conduit de fumée est de 68 °C (voir note de calcul). 

Donc T_iob = 68 °C ≥ T_ig = 53 °C (**Critère 4 : OK**)

