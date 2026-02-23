from flask import Flask, render_template, request, redirect
import pandas as pd
from data_loader import GAF_parser, OBO_parser
from ufficialegraph import GOGraph, _normalise_gaf, _normalise_terms
from analytics import AnalyticsService, CosineSimilarity, JaccardSimilarity
import os
app = Flask(__name__)

# Stato globale vuoto all'avvio
app_state = {'graph': None, 'analytics_engine': None}
'''
@app.route('/')
def index():
    return render_template('index.html')
'''
@app.route('/')
def index():
    # Controlla se il grafo esiste nello stato globale
    is_loaded = app_state['graph'] is not None
    return render_template('index.html', is_loaded=is_loaded)
#------------------------------------------
@app.route('/upload', methods=['POST'])
def upload():
    obo = request.files.get('obo_file')
    gaf = request.files.get('gaf_file')
    
    if not obo or not gaf: 
        return "File mancanti", 400

    obo_path = 'uploaded_ontology.obo'
    gaf_path = 'uploaded_annotations.gaf'
    
    obo.save(obo_path)
    gaf.save(gaf_path)

    try:
        parser_gaf = GAF_parser(gaf_path)
        parser_obo = OBO_parser(obo_path)

        gaf_raw = parser_gaf.load_gaf()
        terms_raw, rels_df = parser_obo.parse_obo()

        app_state['graph'] = GOGraph(_normalise_gaf(gaf_raw), _normalise_terms(terms_raw), rels_df)
        app_state['analytics_engine'] = AnalyticsService(app_state['graph'])

    except Exception as e:
        return f"Errore durante l'elaborazione: {e}", 500

    return redirect('/dashboard')
#----------------------------------------------------
# Rotta dinamica: accetta sia /dashboard che /dashboard/nome_strumento
@app.route('/dashboard')
@app.route('/dashboard/<tool_name>')
def dashboard(tool_name=None):
    if app_state['graph'] is None: 
        return redirect('/')
    # Passa la variabile 'tool_name' all'HTML
    return render_template('dashboard.html', tool_name=tool_name)
#--------------------------------------
@app.route('/search_term', methods=['POST'])
def search_term():
    go_id = request.form.get('go_id', '').strip()
    graph = app_state['graph']
    
    # Recupera i dati grezzi dal parser
    term_raw = graph.query_go_term_detail(go_id, verbose=False)
    
    # 1. RISOLUZIONE CRASH: Interrompi se il termine non esiste
    if not term_raw:
        return f"Termine {go_id} non trovato nel dataset.", 404
    
    t = {}
    # Usa .get() per prevenire ulteriori potenziali KeyError
    t['id'] = term_raw.get('go_id', go_id)
    t['name'] = term_raw.get('name', 'N/A')
    t['definition'] = term_raw.get('definition', 'N/A')
    t['obsoleto'] = term_raw.get('is_obsolete', False)

    # Gestione SINONIMI
    lista_sin = term_raw.get('synonyms', [])
    if lista_sin:
        t['sinonimi'] = ", ".join(lista_sin)
    else:
        t['sinonimi'] = "Nessun sinonimo trovato"

    # 2. RISOLUZIONE BUG VISIVO: Esegui il join della lista sostituti
    lista_sost = term_raw.get('replaced_by', [])
    if lista_sost:
        t['sostituto'] = ", ".join(lista_sost)
    else:
        t['sostituto'] = "Nessun sostituto disponibile"

    # Gestione ANTENATI
    anc_df = graph.query_go_ancestors(go_id)
    if not anc_df.empty:
        t['antenati'] = (anc_df['go_id'] + " - " + anc_df['name']).tolist()
    else:
        t['antenati'] = []

    # Gestione DISCENDENTI
    desc_df = graph.query_go_descendants(go_id)
    if not desc_df.empty:
        t['discendenti'] = (desc_df['go_id'] + " - " + desc_df['name']).tolist()
    else:
        t['discendenti'] = []

    return render_template('term_result.html', term=t)

@app.route('/search_gene', methods=['POST'])
def search_gene():
    if app_state['graph'] is None: return "Errore: Nessun dato caricato. Fai l'upload dei file OBO e GAF dalla Home.", 400

    gene_symbol = request.form.get('gene_symbol', '').strip()
    if not gene_symbol:
        return redirect('/')
        
    graph = app_state['graph']
    
    match = graph.nodes_df[
        (graph.nodes_df['node_type'] == 'gene') & 
        (graph.nodes_df['name'] == gene_symbol)
    ]
    
    if match.empty:
        return render_template('gene_result.html', gene=gene_symbol, annotations=[])
        
    gene_id = match.iloc[0]['node_id']
    annotations = graph.query_gene_annotations(gene_id)
    
    return render_template('gene_result.html', 
                           gene=gene_symbol, 
                           annotations=annotations.to_dict('records'))

@app.route('/similarity', methods=['POST'])
def similarity():
    if app_state['analytics_engine'] is None: return "Errore: Nessun dato caricato. Fai l'upload dei file OBO e GAF dalla Home.", 400

    sim_type = request.form.get('sim_type')
    id1 = request.form.get('id1', '').strip()
    id2 = request.form.get('id2', '').strip()
    metric_choice = request.form.get('metric')
    
    if metric_choice == 'jaccard':
        strategy = JaccardSimilarity()
        metric_name = 'Jaccard'
    else:
        strategy = CosineSimilarity()
        metric_name = 'Cosine'
        
    score = 0.0
    analytics_engine = app_state['analytics_engine']

    try:
        if sim_type == 'term':
            score = analytics_engine.term_similarity(id1, id2, strategy)
        elif sim_type == 'gene':
            score = analytics_engine.gene_similarity(id1, id2, strategy)
    except Exception as e:
        return f"Errore nel calcolo della similarità: {e}", 500
        
    return render_template('similarity_result.html', 
                           sim_type=sim_type, 
                           id1=id1, 
                           id2=id2, 
                           score=score,
                           metric_name=metric_name)

#------------------------------NUOVO STATS
@app.route('/stats', methods=['GET'])
def stats():
    if app_state['graph'] is None or app_state['analytics_engine'] is None: 
        return "Errore: Nessun dato caricato. Fai l'upload dei file OBO e GAF dalla Home.", 400

    graph = app_state['graph']
    analytics_engine = app_state['analytics_engine']
    
    summary_df = analytics_engine.summary()
    top_genes_df = graph.top_annotated_genes(10)
    top_go_df = graph.top_annotated_go_terms(10)
    
    return render_template('stats.html', 
                           summary=summary_df.to_dict('records'),
                           top_genes=top_genes_df.to_dict('records'),
                           top_go=top_go_df.to_dict('records'))

#-----------------------------------------------------neighbourhood
@app.route('/neighbourhood', methods=['POST'])
def neighbourhood():
    if app_state['analytics_engine'] is None: return "Errore: Nessun dato caricato. Fai l'upload dei file OBO e GAF dalla Home.", 400

    go_id = request.form.get('go_id', '').strip()
    depth = request.form.get('depth', 1)
    
    if not go_id:
        return redirect('/')
        
    depth = int(depth)
    neigh_df = app_state['analytics_engine'].neighbourhood_df(go_id, depth)
        
    return render_template('neighbourhood_result.html', 
                           go_id=go_id, 
                           depth=depth,
                           neighbors=neigh_df.to_dict('records'))
#-----------------------------------------------------

if __name__ == '__main__':
    app.run(debug=True, use_reloader=False)