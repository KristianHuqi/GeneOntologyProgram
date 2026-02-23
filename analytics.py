

from abc import ABC, abstractmethod
import numpy as np
import pandas as pd



class OntologyTraversal:
    def __init__(self, relations_df):
        self.rel = relations_df[relations_df["relation_type"] == "is_a"]

    def ancestors(self, go_id):
        visited, frontier = set(), {go_id}
        while frontier:
            parents = set(self.rel[self.rel["child_id"].isin(frontier)]["parent_id"])
            frontier = parents - visited
            visited |= frontier
        visited.add(go_id)
        return visited

    def neighbourhood(self, go_id, depth=1):
        neigh, frontier = {go_id}, {go_id}
        for _ in range(depth):
            parents = set(self.rel[self.rel["child_id"].isin(frontier)]["parent_id"])
            children = set(self.rel[self.rel["parent_id"].isin(frontier)]["child_id"])
            frontier = parents | children
            neigh |= frontier
        return neigh

class SimilarityStrategy(ABC):
    @abstractmethod
    def calculate(self, a, b):
        pass


class CosineSimilarity(SimilarityStrategy):
    def calculate(self, a, b):
        den = float(np.linalg.norm(a) * np.linalg.norm(b))
        return 0.0 if den == 0 else float(a @ b) / den


class JaccardSimilarity(SimilarityStrategy):
    def calculate(self, a, b):
        a, b = a.astype(bool), b.astype(bool)
        inter = np.logical_and(a, b).sum()
        union = np.logical_or(a, b).sum()
        return 0.0 if union == 0 else inter / union


class SimilarityEngine:
    def __init__(self, strategy):
        self.strategy = strategy

    def similarity(self, a, b):
        return self.strategy.calculate(a, b)

class AnalyticsService:
    def __init__(self, graph):
        self.g = graph
        self.tr = OntologyTraversal(graph.relations_df)
    def neighbourhood_df(self, go_id, depth=1):
        ids = self.tr.neighbourhood(go_id, depth)
        return self.g.terms_df[self.g.terms_df["go_id"].isin(ids)][["go_id", "name", "namespace"]]

    def term_similarity(self, t1, t2, strategy):
        a = self.tr.ancestors(t1)
        b = self.tr.ancestors(t2)
        universe = sorted(a | b)
        idx = {x: i for i, x in enumerate(universe)}

        v1 = np.zeros(len(universe), dtype=int)
        v2 = np.zeros(len(universe), dtype=int)

        for x in a:
            v1[idx[x]] = 1
        for x in b:
            v2[idx[x]] = 1

        return SimilarityEngine(strategy).similarity(v1, v2)

    def gene_similarity(self, gene_a, gene_b, strategy):
    
        nodes = self.g.nodes_df
        n_a = nodes[(nodes['node_type'] == 'gene') & (nodes['name'] == gene_a)]
        n_b = nodes[(nodes['node_type'] == 'gene') & (nodes['name'] == gene_b)]
        
        if n_a.empty or n_b.empty:
            raise ValueError("Uno o entrambi i geni non esistono nel database.")
            
        id_a = n_a.iloc[0]['node_id']
        id_b = n_b.iloc[0]['node_id']
        
        
        edges = self.g.edges_df
        
        
        terms_a = set(edges[edges['source'] == id_a]['target']) 
        terms_b = set(edges[edges['source'] == id_b]['target'])
        
        
        tot_annot = sorted(terms_a | terms_b)
        
        
        if not tot_annot:
            return 0.0
            
        idx = {x: i for i, x in enumerate(tot_annot)}
        v1 = np.zeros(len(tot_annot), dtype=int)
        v2 = np.zeros(len(tot_annot), dtype=int)
        
        for x in terms_a: v1[idx[x]] = 1
        for x in terms_b: v2[idx[x]] = 1
            
        return SimilarityEngine(strategy).similarity(v1, v2)
        

    def summary(self):
        nodes, edges = self.g.nodes_df, self.g.edges_df
        rows = [
            ("n_nodes", len(nodes)),
            ("n_edges", len(edges)),
            ("n_go_terms", int((nodes["node_type"] == "GO_term").sum())),
            ("n_genes", int((nodes["node_type"] == "gene").sum())),
            ("n_annotations", int((edges["edge_type"] == "annotates").sum())),
        ]
        return pd.DataFrame(rows, columns=["metric", "value"])