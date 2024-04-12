Pathways
In this chapter we will be going over the pathway file format for MP-BioPath.

OICR Pathway Tab format
This pairwise intearciton tab delimited format reduces the information on the interactions into a computable format.

Simple Example:

5619379    6810158_RLE    1    0
5617644    6810158_RLE    1    0
The file consists of four types of information: entity and interactions. For each interaction there is a parent and a child node and for each pair there is relationship information that specified whether the node is in an "AND" or an "OR" relationship. And if the interaction is a postitive or negative interaction. In this example we have two entities having positibe AND relationship going into a reaction node.

Entities
There is a controlled vocabulary of entity types that can be used as part of a MP-BioPath pathway. These entity types are: gene, rna, dna, protein,If a protein does not have a a preceding RNA entity, one will be created at runtime. And if an RNA entity does not have a preceding genome entity, one will be created at runtime.

A pathway can compose of the following entity types:

Simple Entities
Gene
This is actual gene itself that is used to create RNA

RNA
This is the RNA that is used to create the protein

Protein
protein nodes often interact with other proteins through reactions to effect the functioning of a specific pathway.

Abstract Entities
Set
This entity is used in a pathway when one of many entities could be used to fulfil the same (or similar) function. A set (or Family) is comosed of many "members" or have two sources of the same entity.

Complex
A complex is an entity that is composed of other entites. These entities are known as components of the complex.

Interaction Formats
The first column in contains the parent node. And the second column contains the child node of the interaction. The following are use to specify the information on how the parent child pair interact:

Activating Interaction
1 (Third column)
Inhibition Interaction
-1 (Third column)
AND relationship
0(Fourth column)
Inhibition Interaction
1(Fourth column)
