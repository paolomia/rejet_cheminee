```
Per me si va nella città dolente
per me si va nell'etterno dolore
per me si va tra la perduta gente
Giustizia mosse il mio alto fattore:
dinanzi a me non fuor cose create
se non etterne, e io etterno duro.

Lasciate ogne speranza, voi ch'intrate.
```

J'espère pour toi que tu ne vas pas devoir faire de calculs de rejet de cheminée.

Mais heureusement pour toi, pas mal de travail est déjà fait.

- Télécharge la norme EN 13384-1+A1:2019 (si elle est toujours d'actualité).
- Regarde le `rapport_interne.md` pour avoir un exemple.
- Regarde la `note_calcul/note_de_calcul_AF005215.md` pour avoir un exemple de calcul.

Le conduit de fumée doit respecter plusieurs critères. D'abord lance `calculs.py`, après avoir ajousté les données d'entrée utilisées dans le script.
Le script rédige automatiquement une note de calcul dans `note_de_calcul.md`.

S'il n'y a pas de condensation, tu as de la chance, ton voyage se termine ici.
S'il y a de la condensation, il faut refaire les calculs avec condensation, dans le script `calculs_temperature_cas_condensation.py`. Le script est WIP, donc il faudra probablement ajouster.

Si tu as des hypothèses différentes (apport d'air secondaire, etc.) il faudra regarder directement la norme.

Bon courage