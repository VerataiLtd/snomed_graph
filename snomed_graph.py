import pandas as pd
import networkx as nx
from tqdm.notebook import tqdm
from itertools import groupby
import re
from itertools import pairwise
from typing import Generator, Dict, List, Tuple, Type


class SnomedGraph():

    """
    A class to represent a SNOMED CT release as a Graph, using NetworkX.
    ...

    Attributes
    ----------
    G : nx.DiGraph
        The underlying graph
    """
    

    fsn_typeId = 900000000000003001
    is_a_relationship_typeId = 116680003
    root_concept_id = 138875005

    def __init__(self, G: nx.DiGraph):
        """
        Create a new instance of SnomedGraph from a NetworkX DiGraph object
    
        Args:
            G: A DiGraph created using SnomedGraph.from_rf2() or SnomedGraph.from_serialized().
        Returns:
            self.
        """
        self.G = G
        print(self)

    def __repr__(self):
        return f"SNOMED graph has {self.G.number_of_nodes()} vertices and {self.G.number_of_edges()} edges"

    def __iter__(self):
        return iter(self.G.nodes)       

    def summarise_concept(self, sctid: int) -> None:
        """
        Prints a summary of a concept.
    
        Args:
            sctid: A valid SNOMED Concept ID.
        Returns:
            None
        """        
        parents = self.get_parents(sctid)
        children = self.get_children(sctid)
        inferred_relationships = self.get_inferred_relationships(sctid)
        concept = self.G.nodes[sctid]
        
        print(f"""SCTID:\t\t{sctid}\nFSN:\t\t{concept["fsn"]}\nSynonyms:\t{concept["synonyms"]}""")
        
        print("\nParents:")
        for sctid in parents:
            fsn = self.G.nodes[sctid]["fsn"]
            print(f"""\t{sctid} | {fsn}""")        
            
        print("\nChildren:")        
        for sctid in children:
            fsn = self.G.nodes[sctid]["fsn"]
            print(f"""\t{sctid} | {fsn}""")        
            
        print("\nInferred Relationships:")        
        for group_id, group in inferred_relationships.items():
            print(f"\tGroup {group_id}")
            for g in group:
                print(f"""\t\t---[{g["type"]}]--->\t{g["sctid"]} | {g["fsn"]}""")        
        
    def __get_out_relationships(self, sctid: int) -> Generator[Dict, None, None]:
        for src, tgt in self.G.out_edges(sctid):
            vals = self.G.edges[(src, tgt)]
            node = self.G.nodes[tgt]
            yield {"fsn": node["fsn"], "sctid": tgt, **vals}

    def __get_in_relationships(self, sctid: int) -> Generator[Dict, None, None]:
        for src, tgt in self.G.in_edges(sctid):
            vals = self.G.edges[(src, tgt)]
            node = self.G.nodes[src]
            yield {"fsn": node["fsn"], "sctid": src, **vals}        
    
    def get_concept(self, sctid: int) -> Dict:
        """
        Retrieve all attributes for a given concept.
    
        Args:
            sctid: A valid SNOMED Concept ID.
        Returns:
            A dictionary containing all the attributes of the concept.
        """
        return {"sctid": sctid, **self.G.nodes[sctid]}

    def get_parents(self, sctid: int) -> List[int]:
        """Retrieve all parents of a given concept.
    
        Args:
            sctid: A valid SNOMED Concept ID.
        Returns:
            A list containing the SCTIDs of all parents.
        """        
        return [
            p["sctid"]
            for p in self.__get_out_relationships(sctid) 
            if p["group"] == 0 and p["type_id"] == SnomedGraph.is_a_relationship_typeId
        ]

    def get_children(self, sctid: int) -> List[int]:
        """
        Retrieve all children of a given concept.
    
        Args:
            sctid: A valid SNOMED Concept ID.
        Returns:
            A list containing the SCTIDs of all children.
        """
        return [
            c["sctid"]
            for c in self.__get_in_relationships(sctid) 
            if c["group"] == 0 and c["type_id"] == SnomedGraph.is_a_relationship_typeId
        ]

    def get_descendants(self, sctid: int, steps_removed: int = None) -> List[int]:
        """
        Retrieve descendants of a given concept.
    
        Args:
            sctid: A valid SNOMED Concept ID.
            steps_removed: The number of levels down in the hierarchy to go.  
                           (1 => children; 2 => children + grandchildren, etc)
                           if None then all children are retrieved.
        Returns:
            A list containing the SCTIDs of all descendants.
        """
        if steps_removed is None:
            steps_removed = 99999
        elif steps_removed <= 0:
            raise AssertionError("steps_removed must be > 0 or None")
        children = self.get_children(sctid)
        descendants = set(children)
        if steps_removed > 1:
            for c in children:
                descendants = descendants.union(
                    [d for d in self.get_descendants(c, steps_removed-1)]
                )
        return descendants

    def get_ancestors(self, sctid: int, steps_removed: int = None) -> List[int]:
        """
        Retrieve ancestors of a given concept.
    
        Args:
            sctid: A valid SNOMED Concept ID.
            steps_removed: The number of levels up in the hierarchy to go.  
                           (1 => parents; 2 => parents + grandparents, etc)
                           if None then all parents are retrieved.
        Returns:
            A list containing the SCTIDs of all descendants.
        """        
        if steps_removed is None:
            steps_removed = 99999
        elif steps_removed <= 0:
            raise AssertionError("steps_removed must be > 0 or None")
        parents = self.get_parents(sctid)
        ancestors = set(parents)
        if steps_removed > 1:
            for p in parents:
                ancestors = ancestors.union(
                    [a for a in self.get_ancestors(p, steps_removed-1)]
                )
        return ancestors.difference({SnomedGraph.root_concept_id})

    def get_neighbourhood(self, sctid: int, steps_removed: int = 1) -> List[int]:
        """
        Retrieve neighbours of a given concept.
        Neighbours include ancestors, descendants and cousins up to the given degree.
    
        Args:
            sctid: A valid SNOMED Concept ID.
            steps_removed: The number of steps up or down in the hierarchy to go.
                           Defaults to 1 (parents + children).
        Returns:
            A list containing the SCTIDs of all neighbours.
        """            
        assert steps_removed > 0
        parents = self.get_parents(sctid)
        children = self.get_children(sctid)
        neighbourhood = set(parents).union(children)
        if steps_removed > 1:
            for n in list(neighbourhood):
                neighbourhood = neighbourhood.union(
                    [p for p in self.get_neighbourhood(n, steps_removed-1)]
                )
                neighbourhood = neighbourhood.union(
                    [c for c in self.get_neighbourhood(n, steps_removed-1)]
                )
            neighbourhood = neighbourhood.difference([sctid])
        return neighbourhood.difference({SnomedGraph.root_concept_id})
    
    def get_inferred_relationships(self, sctid: int) -> Dict:
        """
        Retrieve the "Inferred Relationships" for a concept.
        Inferred relationships exist in groups.  Together with the parents of a concept 
        the set of inferred relationship groups constitute a unique specification of the
        concept.
    
        Args:
            sctid: A valid SNOMED Concept ID.
        Returns:
            A dictionary with keys = group IDs and values = list of inferred relationships
        """          
        key_ = lambda r: r["group"]
        it = groupby(sorted(self.__get_out_relationships(sctid), key=key_), key=key_)
        return {
            group: [
                {"type": r["type"], "fsn": r["fsn"], "sctid": r["sctid"]} 
                for r in relationships
                if r["type_id"] != SnomedGraph.is_a_relationship_typeId
            ] 
            for group, relationships in it
        }

    def find_path(self, sctid1: int, sctid2: int, print_: bool = False) -> List[Tuple]:
        """
        Returns details of any path that exists between two concepts.
        The path considers all relationship types but limits the results to true ancestors
        or descentants - i.e. concepts that are "cousins" of one another will not result in
        a returned path.
    
        Args:
            sctid1: A valid SNOMED Concept ID.
            sctid2: A valid SNOMED Concept ID.
            print_: Whether to print the full relationship as a string.
        Returns:
            A list of tuples of the form (source, relationship_type, target).  The first
            'source' element of the first tuple will be one of sctid1 or sctid2.  The last
            'target' element of the last tuple will be the other sctid1 or sctid2.
        """         
        path = []
        if nx.has_path(self.G, sctid1, sctid2):
            nodes = nx.shortest_path(self.G, sctid1, sctid2)
        elif nx.has_path(self.G, sctid2, sctid1):
            nodes = nx.shortest_path(self.G, sctid2, sctid1)
        else:
            nodes = []
        for src, tgt in pairwise(nodes):
            r = self.G.edges[(src, tgt)]
            path.append((src, r["type"], tgt))
        if print_:
            src = path[0][0]
            str_ = f"({self.get_concept(src)['fsn']})"
            for _, rel, tgt in path:
                str_ += f"---[{rel}]--->"
                str_ += f"({self.get_concept(tgt)['fsn']})"
            print(str_)
        return path

    def save(self, path: str) -> None:
        """
        Save this SnomedGraph
    
        Args:
            path: path + name of the file to save to.
        Returns:
            None
        """          
        nx.write_gml(self.G, path)

    @staticmethod
    def from_serialized(path: str):
        """
        Load a SnomedGraph from a serialization.
    
        Args:
            path: path + name of the file to save to.
        Returns:
            A SnomedGraph
        """            
        G = nx.read_gml(path, destringizer=int)
        return SnomedGraph(G)    
        
    @staticmethod
    def from_rf2(path: str):
        """
        Create a SnomedGraph from a SNOMED RF2 release path.
    
        Args:
            path: Path to RF2 release folder.
        Returns:
            A SnomedGraph
        """          
        if path[-1] == "/":
            path = path[:-1]
        release_date_pattern = r'\d{8}'
        match = re.search(release_date_pattern, path)
        try:                        
            release_date = match.group(0)
        except AttributeError:
            raise AssertionError(
                f"The path does not appear to contain a valid SNOMED CT Release Format name."
            )
        else:
            # Load relationships
            relationships_df = pd.read_csv(
                f"{path}/Snapshot/Terminology/sct2_Relationship_Snapshot_INT_{release_date}.txt", 
                delimiter="\t"
            )
            relationships_df = relationships_df[relationships_df.active == 1]
            
            # Load concepts
            concepts_df = pd.read_csv(
                f"{path}/Snapshot/Terminology/sct2_Description_Snapshot-en_INT_{release_date}.txt", 
                delimiter="\t"
            )
            concepts_df = concepts_df[concepts_df.active == 1]
            concepts_df.set_index("conceptId", inplace=True)

            # Create relationships type lookup
            relationship_types = concepts_df.loc[relationships_df.typeId.unique()]
            relationship_types = relationship_types[relationship_types.typeId == SnomedGraph.fsn_typeId]
            relationship_types = relationship_types.term.to_dict()

            # Initialise the graph
            n_concepts = concepts_df.shape[0]
            n_relationships = relationships_df.shape[0]
            print(f"{n_concepts} terms and {n_relationships} relationships were found in the release.")
            G = nx.DiGraph()            

            # Create relationships
            print("Creating Relationships...")
            for r in tqdm(relationships_df.to_dict(orient="records")):
                G.add_edge(        
                    r["sourceId"], 
                    r["destinationId"], 
                    group=r["relationshipGroup"], 
                    type=relationship_types[r["typeId"]],
                    type_id=r["typeId"]
                )

            # Add concepts            
            print("Adding Concepts...")
            for sctid, rows in tqdm(concepts_df.groupby(concepts_df.index)):
                synonyms = [row.term for _, row in rows.iterrows() if row.typeId != SnomedGraph.fsn_typeId]
                try:
                    fsn = rows[rows.typeId == SnomedGraph.fsn_typeId].term.values[0]
                except IndexError:
                    fsn = synonyms[0]
                    synonyms = synonyms[1:]
                    print(f"Concept with SCTID {sctid} has no FSN. Using synonym '{fsn}' instead.")
                G.add_node(sctid, fsn=fsn, synonyms=synonyms)

            # Remove isolates
            G.remove_nodes_from(list(nx.isolates(G)))
            
            # Initialise class            
            return SnomedGraph(G)
            