
import pandas as pd
from data_loader import GAF_parser, OBO_parser
import numpy as np


_GAF_COLUMN_MAP = {
    'DB object id':         'DB_object_id',
    'DB object symbol':     'DB_object_symbol',
    'GO ID':                'GO_ID',
    'Evidence code':        'Evidence_code',
    'DB object name':       'DB_object_name',
    'DB object synonym':    'DB_object_synonym',
    'DB object type':       'DB_object_type',
    'Assigned by':          'Assigned_by',
    'Annotation Extension': 'Annotation_extension',
    'Gene product form id': 'Gene_product_form_id',
    'With (or) From':       'With_from',
}

_TERMS_COLUMN_MAP = {
    'def': 'definition',
}


def _normalise_gaf(df):
    return df.rename(columns=_GAF_COLUMN_MAP)


def _normalise_terms(df):
    df = df.rename(columns=_TERMS_COLUMN_MAP)
    for col, default in [
        ('namespace',   ''),
        ('definition',  ''),
        ('is_obsolete', False),
        ('replaced_by', ''),
        ('consider',    ''),
    ]:
        if col not in df.columns:
            df[col] = default
    return df


class GraphBuilder:

    def __init__(self, gaf_df, terms_df, relations_df,
                 include_obsolete=False, evidence_codes=None, aspects=None):
        self.gaf_df           = gaf_df
        self.terms_df         = terms_df
        self.relations_df     = relations_df
        self.include_obsolete = include_obsolete
        self.evidence_codes   = evidence_codes
        self.aspects          = aspects
        self.nodes_df         = pd.DataFrame()
        self.edges_df         = pd.DataFrame()
        self._build()

    def _build(self):
        terms = self.terms_df.copy()
        if not self.include_obsolete:
            terms = terms[~terms['is_obsolete']].copy()

        go_nodes = terms[['go_id', 'name', 'namespace', 'is_obsolete']].copy()
        go_nodes.columns = ['node_id', 'name', 'namespace', 'is_obsolete']
        go_nodes['node_type'] = 'GO_term'

        gaf_filtered = self.gaf_df.copy()
        if self.evidence_codes:
            gaf_filtered = gaf_filtered[gaf_filtered['Evidence_code'].isin(self.evidence_codes)]
        if self.aspects:
            gaf_filtered = gaf_filtered[gaf_filtered['Aspect'].isin(self.aspects)]

        gene_nodes = (
            gaf_filtered[['DB_object_id', 'DB_object_symbol', 'DB']]
            .drop_duplicates(subset='DB_object_id')
            .rename(columns={
                'DB_object_id':     'node_id',
                'DB_object_symbol': 'name',
                'DB':               'namespace',
            })
        )
        gene_nodes['node_type']   = 'gene'
        gene_nodes['is_obsolete'] = False

        self.nodes_df = (
            pd.concat([go_nodes, gene_nodes], ignore_index=True)
            .drop_duplicates(subset='node_id')
            .reset_index(drop=True)
        )

        valid_ids = set(self.nodes_df['node_id'])

        onto_edges = (
            self.relations_df[
                self.relations_df['child_id'].isin(valid_ids) &
                self.relations_df['parent_id'].isin(valid_ids)
            ]
            .copy()
            .rename(columns={
                'child_id':      'source',
                'parent_id':     'target',
                'relation_type': 'edge_type',
            })
        )
        onto_edges['evidence'] = ''
        onto_edges['aspect']   = ''

        ann_edges = gaf_filtered[['DB_object_id', 'GO_ID', 'Evidence_code', 'Aspect']].copy()
        ann_edges.columns = ['source', 'target', 'evidence', 'aspect']
        ann_edges = ann_edges[ann_edges['target'].isin(valid_ids)]
        ann_edges['edge_type'] = 'annotates'

        self.edges_df = pd.concat([onto_edges, ann_edges], ignore_index=True)
        print(f'[GRAPH] nodes: {len(self.nodes_df):,}  |  edges: {len(self.edges_df):,}')


class GraphQuery(GraphBuilder):

    def query_gene_annotations(self, gene_id):
        edges = self.edges_df[
            (self.edges_df['source'] == gene_id) &
            (self.edges_df['edge_type'] == 'annotates')
        ]
        result = edges.merge(
            self.nodes_df[['node_id', 'name', 'namespace']],
            left_on='target', right_on='node_id', how='left',
        )[['target', 'name', 'namespace', 'evidence', 'aspect']]
        result.columns = ['GO_ID', 'GO_name', 'namespace', 'evidence', 'aspect']
        return result.reset_index(drop=True)

    def query_go_genes(self, go_id):
        edges = self.edges_df[
            (self.edges_df['target'] == go_id) &
            (self.edges_df['edge_type'] == 'annotates')
        ]
        result = edges.merge(
            self.nodes_df[['node_id', 'name']],
            left_on='source', right_on='node_id', how='left',
        )[['source', 'name', 'evidence']]
        result.columns = ['gene_id', 'gene_name', 'evidence']
        return result.reset_index(drop=True)

    def query_go_ancestors(self, go_id, max_depth=10):
        visited, frontier = set(), {go_id}
        for _ in range(max_depth):
            if not frontier:
                break
            parents = self.relations_df[self.relations_df['child_id'].isin(frontier)]['parent_id']
            new = set(parents) - visited - frontier
            visited |= frontier
            frontier = new
        return (
            self.terms_df[self.terms_df['go_id'].isin(visited - {go_id})]
            [['go_id', 'name', 'namespace']]
            .reset_index(drop=True)
        )

    def query_go_descendants(self, go_id, max_depth=10):
        visited, frontier = set(), {go_id}
        for _ in range(max_depth):
            if not frontier:
                break
            children = self.relations_df[self.relations_df['parent_id'].isin(frontier)]['child_id']
            new = set(children) - visited - frontier
            visited |= frontier
            frontier = new
        return (
            self.terms_df[self.terms_df['go_id'].isin(visited - {go_id})]
            [['go_id', 'name', 'namespace']]
            .reset_index(drop=True)
        )

    def query_go_term_detail(self, go_id, verbose=True):
        row = self.terms_df[self.terms_df['go_id'] == go_id]
        if row.empty:
            print(f'[WARN] {go_id} not found.')
            return {}
        row = row.iloc[0]

        is_a_parents = self.relations_df[
            (self.relations_df['child_id'] == go_id) &
            (self.relations_df['relation_type'] == 'is_a')
        ]['parent_id'].tolist()

        other_rels = self.relations_df[
            (self.relations_df['child_id'] == go_id) &
            (self.relations_df['relation_type'] != 'is_a')
        ][['relation_type', 'parent_id']].values.tolist()

        detail = {
            'go_id':         row['go_id'],
            'name':          row.get('name', ''),
            'namespace':     row.get('namespace', ''),
            'definition':    row.get('definition', ''),
            'synonyms':      row.get('synonyms', []),
            'is_obsolete':   row.get('is_obsolete', False),
            'is_a':          is_a_parents,
            'relationships': [{'type': r[0], 'target': r[1]} for r in other_rels],
            'replaced_by':   [c for c in str(row.get('replaced_by', '')).split('; ') if c],
            'consider':      [c for c in str(row.get('consider', '')).split('; ') if c],
        }
        if verbose:
            self._print_term_detail(detail)
        return detail

    @staticmethod
    def _print_term_detail(d):
        print('─' * 60)
        print(f"  id        : {d['go_id']}")
        print(f"  name      : {d['name']}")
        print(f"  namespace : {d['namespace']}")
        print(f"  obsolete  : {d['is_obsolete']}")
        print(f"  def       : {d['definition']}")
        if d['is_a']:
            print('  is_a      :')
            for pid in d['is_a']:
                print(f'    {pid}')
        else:
            print('  is_a      : (none)')
        if d['relationships']:
            print('  relations :')
            for r in d['relationships']:
                print(f"    {r['type']:20s}  {r['target']}")
        if d['replaced_by']:
            print('  replaced_by:')
            for rid in d['replaced_by']:
                print(f'    {rid}')
        print()


class GraphStats(GraphQuery):

    def get_stats(self):
        go_nodes    = (self.nodes_df['node_type'] == 'GO_term').sum()
        gene_nodes  = (self.nodes_df['node_type'] == 'gene').sum()
        ns_counts   = self.nodes_df[self.nodes_df['node_type'] == 'GO_term']['namespace'].value_counts()
        edge_counts = self.edges_df['edge_type'].value_counts()
        rows = [
            ('total_nodes', len(self.nodes_df)),
            ('total_edges', len(self.edges_df)),
            ('go_terms',    int(go_nodes)),
            ('genes',       int(gene_nodes)),
        ]
        for ns, cnt in ns_counts.items():
            rows.append((f'namespace:{ns}', int(cnt)))
        for et, cnt in edge_counts.items():
            rows.append((f'edges:{et}', int(cnt)))
        return pd.DataFrame(rows, columns=['metric', 'value'])

    def top_annotated_genes(self, n=10):
        ann    = self.edges_df[self.edges_df['edge_type'] == 'annotates']
        counts = ann.groupby('source').size().reset_index(name='n_annotations')
        counts = counts.merge(
            self.nodes_df[['node_id', 'name']], left_on='source', right_on='node_id', how='left',
        )[['source', 'name', 'n_annotations']]
        counts.columns = ['gene_id', 'gene_name', 'n_annotations']
        return counts.sort_values('n_annotations', ascending=False).head(n).reset_index(drop=True)

    def top_annotated_go_terms(self, n=10):
        ann    = self.edges_df[self.edges_df['edge_type'] == 'annotates']
        counts = ann.groupby('target').size().reset_index(name='n_genes')
        counts = counts.merge(
            self.nodes_df[['node_id', 'name', 'namespace']], left_on='target', right_on='node_id', how='left',
        )[['target', 'name', 'namespace', 'n_genes']]
        counts.columns = ['GO_ID', 'GO_name', 'namespace', 'n_genes']
        return counts.sort_values('n_genes', ascending=False).head(n).reset_index(drop=True)


class GOGraph(GraphStats):
    
    def go_adjacency_matrix(self):
        rel = self.relations_df[self.relations_df["relation_type"] == "is_a"]
        terms = sorted(set(rel["child_id"]) | set(rel["parent_id"]))
        A = np.zeros((len(terms), len(terms)), dtype=int)

        for i, child in enumerate(terms):
            parents = rel[rel["child_id"] == child]["parent_id"]
            for parent in parents:
                j = terms.index(parent)
                A[i, j] = 1

        return A, terms

    def gene_go_matrix(self):
        genes = sorted(self.nodes_df[self.nodes_df["node_type"] == "gene"]["node_id"])
        terms = sorted(self.nodes_df[self.nodes_df["node_type"] == "GO_term"]["node_id"])

        M = np.zeros((len(genes), len(terms)), dtype=int)

        ann = self.edges_df[self.edges_df["edge_type"] == "annotates"]

        for _, row in ann.iterrows():
            if row["source"] in genes and row["target"] in terms:
                i = genes.index(row["source"])
                j = terms.index(row["target"])
                M[i, j] = 1

        return M, genes, terms

    
if __name__ == '__main__':

    FILEPATH_GAF = 'data/goa_human.gaf'
    FILEPATH_OBO = 'data/go-basic.obo'

    parser_gaf = GAF_parser(FILEPATH_GAF)
    parser_obo = OBO_parser(FILEPATH_OBO)

    gaf_raw            = parser_gaf.load_gaf()
    terms_raw, rels_df = parser_obo.parse_obo()
    

    gaf_df   = _normalise_gaf(gaf_raw)
    terms_df = _normalise_terms(terms_raw)

    graph = GOGraph(
        gaf_df       = gaf_df,
        terms_df     = terms_df,
        relations_df = rels_df,
    )

    
    print('\n── Stats ──')
    print(graph.get_stats().to_string(index=False))

    go = 'GO:0003723'
    print(f'\n── Detail for {go} ──')
    graph.query_go_term_detail(go)

    gene = 'A0A024RBG1'
    print(f'\n── GO annotations for {gene} ──')
    print(graph.query_gene_annotations(gene).to_string(index=False))

    print(f'\n── Ancestors of {go} ──')
    print(graph.query_go_ancestors(go).to_string(index=False))

    print(f'\n── Descendants of {go} ──')
    print(graph.query_go_descendants(go).to_string(index=False))

    print('\n── Top 10 most annotated genes ──')
    print(graph.top_annotated_genes().to_string(index=False))

    print('\n── Top 10 most populated GO terms ──')
    print(graph.top_annotated_go_terms().to_string(index=False))
