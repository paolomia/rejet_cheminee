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

Les conduits sont en tôle métallique non isolé, 1 mm d'épaisseur. Le conduit situé à l'intérieur du batimen est dans un environnement chauffé (15 °C).

Un volet de régulation est installé dans le conduit de fumée.

L'adresse de l'installation est Havighorster Weg 8C, 21031 Hamburg, Germania. Le site est situé à plus de 20 km de la mer et à une altitude de 41 m au-dessus du niveau de la mer.

D'après tests interns, nous pouvons estimer un débit massique pour la fumée d'environs 0.009 kg/s et une température de fumée à la sortie du foyer de 300 °C.

Je n'ai pas de données concernant la pression maximale supportée par le conduit de fumée.

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

# 2. Calculs

## 2.1 Pression - Critères 1-2
Nous avons les plage de pression de fonctionnement du bruleur, mais pas de l'ensemble bruleur + foyer. Cependant, j'estime que nous pouvons quand même justifier l'adhérence du conduit aux critères de pression imposés.

Le conduit de fumée dispose d'un volet réglable pour le bon fonctionnement du bruleur. Soit **P_R_volet** la perte de charge en Pa engendrée par le volet.

D'après les calculs, nous obtenons les valeurs suivant (estimation un peu grossière, pour presenter l'argument) :

P_ZO ≈ -8 Pa + P_R_volet
P_ZOmin ≈ -43 Pa + P_R_volet

Essentiellement, étant donné que la vitesse de la fumée est extrêmement faible, la perte de charge engendrée par le conduit est négligeable, avec une valeur de tirage bien plus importante (40 Pa) qui peut être compensée par le volet.

Nous avons le graphique suivant pour les plages de fonctionnement du bruleur :
![img.png](img.png)

donc entre -50 Pa et 300 Pa.

Il faut trouver comment le presenter auprès du client, mais à mon avis, avec ces données, nous respectons largement les critères de pression 1-2.

Nous avons une cheminée avec une perte de charge adjustable, à partir de -8 Pa avec une plage de pression qui varie de moins de 40 Pa

Le bruleur supporte une plage de variation bien plus importante. Si les critères de pression n'étaient pas ici respectés, alors ça ne serait pas la cheminée à devoir être remise en question, mais le système bruleur - foyer échangeur. 

## 2.2 Pression - Critère 3

Je n'ai pas d'élément concernant la pression maximale supportée par nos conduits de fumée. 
Réalistiquement, j'imagine que nos conduits peuvent supporter une pression d'au moins 300 Pa, ce qui est la pression maximale de fonctionnement du bruleur ?

## 2.3 Température - Critère 4

La température de rosée est estimée à 54 °C. Même si la température de la fumée à la sortie est autour de ~110 °C, la condution thermique avec l'air ambiant est bien plus importante, ce qui implique que nous avons sans doute de la condensation dans le conduit.

Cela nous oblige à réaliser le calcul avec condensation (Section 8, pg 82) ce qui est un calcul encore plus tédieux, qui demande de découper le conduit en segment, et pour chacun resoudre un système d'équations qui décrivent le processus de transfert thermique du à la condensation.
Il faut que la température de la paroi intérieure à la sortie du conduit de fumée soit ≥ 0 °C.

Nous sommes vraiment à la limite, en variant légèrement les hypothèses on peut tomber legèrement en dessous ou au dessus de la température de 0 °C.

