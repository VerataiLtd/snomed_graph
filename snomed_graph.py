import pandas as pd
import networkx as nx
from tqdm.notebook import tqdm
from itertools import groupby
import re


class SnomedGraph():

    fsn_typeId = 900000000000003001
    synonym_typeId = 900000000000013009

    def __init__(self, G):        
        assert isinstance(G, nx.DiGraph)
        self.G = G
        print(self)

    def __repr__(self):
        return f"SNOMED graph has {self.G.number_of_nodes()} vertices and {self.G.number_of_edges()} edges"

    def __iter__(self):
        return iter(self.G.nodes)       

    def summarise_concept(self, sctid):
        parents = self.get_parents(sctid)
        children = self.get_children(sctid)
        inferred_relationships = self.get_inferred_relationships(sctid)
        concept = self.G.nodes[sctid]
        
        print(f"""SCTID:\t\t{sctid}\nFSN:\t\t{concept["name"]}\nSynonyms:\t{concept["synonyms"]}""")
        
        print("\nParents:")
        for p in parents:
            print(f"""\t{p["concept_id"]} | {p["name"]}""")        
            
        print("\nChildren:")        
        for c in children:
            print(f"""\t{c["concept_id"]} | {c["name"]}""")        
            
        print("\nInferred Relationships:")        
        for group_id, group in inferred_relationships.items():
            print(f"\tGroup {group_id}")
            for g in group:
                print(f"""\t\t---[{g["type"]}]--->\t{g["concept_id"]} | {g["name"]}""")        
        
    def __get_out_relationships(self, sctid):
        for src, tgt in self.G.out_edges(sctid):
            vals = self.G.edges[(src, tgt)]
            node = self.G.nodes[tgt]
            yield {"name": node["name"], "concept_id": tgt, **vals}

    def __get_in_relationships(self, sctid):
        for src, tgt in self.G.in_edges(sctid):
            vals = self.G.edges[(src, tgt)]
            node = self.G.nodes[src]
            yield {"name": node["name"], "concept_id": src, **vals}        
    
    def get_concept(self, sctid):
        return {"concept_id": sctid, **self.G.nodes[sctid]}

    def get_parents(self, sctid):
        return [
            {"name": p["name"], "concept_id": p["concept_id"]}
            for p in self.__get_out_relationships(sctid) 
            if p["group"] == 0
        ]

    def get_children(self, sctid):
        return [
            {"name": c["name"], "concept_id": c["concept_id"]}
            for c in self.__get_in_relationships(sctid) 
            if c["group"] == 0
        ]

    def get_inferred_relationships(self, sctid):
        key_ = lambda r: r["group"]
        it = groupby(sorted(self.__get_out_relationships(sctid), key=key_), key=key_)
        return {
            group: [
                {"type": r["type"], "name": r["name"], "concept_id": r["concept_id"]} 
                for r in relationships
            ] 
            for group, relationships in it
            if group > 0
        }        

    def save(self, path):
        nx.write_gml(G, path)

    @staticmethod
    def from_serialized(path):
        G = nx.read_gml(path, destringizer=int)
        return SnomedGraph(G)    
        
    @staticmethod
    def from_rf2(path):
        if path[-1] == "/":
            path = path[:-1]
        release_date_pattern = r'\d{8}'
        match = re.search(release_date_pattern, path)
        try:                        
            release_date = match.group(0)
        except AttributeError:
            raise AssertionError(f"The path does not appear to contain a valid SNOMED CT Release Format name.")        
        else:
            # Load relationships
            relationships_df = pd.read_csv(f"{path}/Snapshot/Terminology/sct2_Relationship_Snapshot_INT_{release_date}.txt", delimiter="\t")
            relationships_df = relationships_df[relationships_df.active == 1]
            
            # Load concepts
            concepts_df = pd.read_csv(f"{path}/Snapshot/Terminology/sct2_Description_Snapshot-en_INT_{release_date}.txt", delimiter="\t")
            concepts_df = concepts_df[concepts_df.active == 1]
            concepts_df.set_index("conceptId", inplace=True)

            # Create relationships type lookup
            relationship_types = concepts_df.loc[relationships_df.typeId.unique()]
            relationship_types = relationship_types[relationship_types.typeId == fsn_typeId]
            relationship_types = relationship_types.term.to_dict()

            print(f"{concepts_df.shape[0]} terms and {relationships_df.shape[0]} relationships were found in the release.")
            G = nx.DiGraph()            

            # Create relationships
            print("Creating Relationships...")
            for r in tqdm(relationships_df.to_dict(orient="records")):
                G.add_edge(        
                    r["sourceId"], 
                    r["destinationId"], 
                    group=r["relationshipGroup"], 
                    type=relationship_types[r["typeId"]]
                )

            # Add concepts            
            print("Adding Concepts...")
            for concept_id, rows in tqdm(concepts_df.groupby(concepts_df.index)):
                synonyms = [row.term for _, row in rows.iterrows() if row.typeId != fsn_typeId]
                try:
                    fsn = rows[rows.typeId == fsn_typeId].term.values[0]
                except IndexError:
                    fsn = synonyms[0]
                G.add_node(concept_id, name=fsn, synonyms=synonyms)

            # Initialise class            
            return SnomedGraph(G)
            