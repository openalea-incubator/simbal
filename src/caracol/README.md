# Code organisation

## Parameters

- Create a file to read/write the parameters needed for the simulation.
With default values and unit (see walter?).

- Separate the definition of parameters from the code.
    * in a parameter file (resp. openalea)
    * in a test, explicitly in python
    * in visualea, in a given workflow

### Basic
- Set the parameters directly in a test.

## Architecture

Define global simulation programs of several plant architectures (various varieties) (a list) with an above and below ground part.
Each part is represented as an MTG. But may have different scale and properties for the above and below ground part.
As well as for between plants.

## Allocation CN

Manage the allocation (shoot & root) of C and N.

## Growth

Define the growth of the plants considering C and N allocation.

## Geometry

Reconstruct in 3D the plants' geometry 

## Environment

Compute light interception and soil N content

## Senescence

### Algorithm

**TDOD** : translate in english

Internal function :

- function: *leaf_status* 
Détermine l'état de chaque feuille du MTG et je leurs attribue des propriétés ;
Si teneur en C insuffisante face à la demande --> dC négatif
alors growth = false
    storage = false
    sink = true
    source = false
Sinon la teneur en C est suffisante face à la demande --> dC positif
alors growth = true
    storage = true
    sink = false
    source = true

- C_offer (photosynthesis)
Détermine la participation des feuilles 'source' du MTG au pool C et je met à jour leur propriété teneur en C ;
Plus la lumière interceptée (PAR) de la feuille est élevé, plus son poids est important --> C_leaf_export
C_leaf_actual_content = C__potential_leaf_offer - C_leaf_export

- Je fais le bilan du pool C ;
Soit la somme du carbone exporté de l'ensemble des feuilles du MTG --> C_pool_offer = sum(C_leaf_export)
Si C_pool_offer > 0
alors empty = false, il existe une réserve en carbone disponible
Sinon empty = true, le pool est vide

- Je vérifie que la réserve en C est suffisante face à la demande ;
dA = dC + C_pool_offer
Si le pool est insuffisant --> dA négatif
alors stress = true
    growth = false
    senescence_initiation = true
    storage = 0
    C_actual_deficit = dA
Sinon le pool est suffisant --> dA positif
alors stress = false
    growth = true
    senescence_initiation = false
    storage = dA
    Je sélectionne les feuilles 'puits' du MTG et je met à jour leur propriété --> C_actual_uptake = C_potential_leaf_demand

- Je développe le mécanisme de sénescence (leaf_senescence function) ;
Je met à jour les propriétés de l'ensemble des feuilles du MTG :
J'évalue la teneur en N --> N_leaf_actual_content = (N_potential_plant_offer / plante_actual_surface) * leaf_area

Je sélectionne les feuilles 'source' du MTG :
Je sélectionne la feuille 'source' avec l'interception de la lumière (PAR) la plus faible
Je met à jour ses propriété et annonce sa sénescence et les seuils d'initiation de la sénescence --> 'source' = 'sénescence', N_senescence_threshold = N_leaf_actual_content, C_senescence_threshold = C_leaf_actual_content

Je détermine la participation des feuilles 'sénescence' du MTG au pool N et je met à jour leur propriété teneur en N
Plus la lumière interceptée (PAR) de la feuille est élevé, plus son poids est important --> N_leaf_remobilized = w * N_leaf_actual_content
N_leaf_senescent_content = N_leaf_actual_content - N_leaf_remobilized 

Je détermine la valeur du pool N :
Je totalise l'ensemble des feuilles 'sénescence' du MTG --> sum(N_leaf_remobilized)

Je détermine la demande en N de la plante :
Je sélectionne les feuilles 'source' et 'puit' du MTG
J'évalue la demande en azote de ces feuilles et leur attribue cette propriété --> N_leaf_deficit = N_leaf_actual_content - leaf_N_potential_offer
Soit la somme de l'azote demandé par les feuilles 'source' et 'puit' du MTG --> N_plant_demand = sum(N_leaf_deficit)

Je vérifie que la réserve est suffisante face à la demande en N :
N_plant_hunger = N_pool_offer - N_plant_demand
Si pool suffisant --> N_plant_hunger positif
alors N_plant_remobilized = N_plant_demand
    N_plant_available = N_plant_hunger
    La demande est entièrement satisfaite et je met à jour leur propriété teneur en N --> N_leaf_actual_content = leaf_N_potential_offer
Sinon pool insuffisant --> N_plant_hunger négatif
    Je détermine l'importance de chaque feuille 'puit' et 'source' en azote importé selon sa lumière interceptée (PAR) et je met à jour leur propriété teneur en N
    Plus la lumière interceptée (PAR) de la feuille est élevé, plus son poids est important --> N_leaf_import = w * N_pool_offer
    N_leaf_deficit = N_leaf_actual_content + N_leaf_import
    J'annonce la chute de l'ensemble des feuilles 'sénescence' et je met à jour leur teneur en N  --> N_leaf_fall = N_leaf_actual_content - N_leaf_remobilized
    Leaf_fall = true
    


Main function :

- Je détermine l'état de croissance de la feuille
Si pas en croissance --> growth = false
alors je vérifie s'il y a une réserve en C
    Si un pool est disponible --> empty = false
    alors je vérifie que ce pool est suffisant
        Si pool suffisant --> stress = false
        Je met à jour les propriétés de la feuille --> C_actual_uptake = C_potential_leaf_demand
    Sinon pool insuffisant --> stress = true
    J'annonce la sénescence d'une feuille 'source' du MTG --> leaf_senescence()
Sinon croissance optimale --> growth = true
Je met à jour les propriétés de l'ensemble des feuilles 'puit' et 'source' du MTG --> C_actual_uptake = C_potential_leaf_demand














