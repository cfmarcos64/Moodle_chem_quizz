# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 12:45:52 2026

@author: Carlos Fernandez
"""

"""
Interactive Molecular Editor Module (Canvas / streamlit-ketcher)
Allows drawing chemical structures and automatically extracting SMILES, formulas, 
PubChem IUPAC names, and generating Moodle questions in a two-column layout.
"""

import streamlit as st
import io
import base64
import urllib.request
import json

# --- Dependency Verification ---
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, FindMolChiralCenters
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

try:
    from streamlit_ketcher import st_ketcher
    KETCHER_AVAILABLE = True
except ImportError:
    KETCHER_AVAILABLE = False

# Bilingual interface texts
TEXTS_EDITOR = {
    "es": {
        "title": "Editor Molecular Interactivo (Canvas) ✍️",
        "desc": "Dibuja estructuras químicas en el lienzo, pulsa 'Apply' en el editor y obtén al instante su análisis, nombre IUPAC y exportación a Moodle.",
        "canvas_label": "Lienzo de Dibujo Químico",
        "no_ketcher": "La librería `streamlit-ketcher` no está instalada.",
        "no_rdkit": "RDKit no está disponible en este entorno.",
        "empty_canvas": "Dibuja una molécula en el lienzo superior y pulsa 'Apply' para ver su análisis.",
        "molecular_formula": "Fórmula Molecular",
        "molecular_weight": "Peso Molecular",
        "canonical_smiles": "SMILES Canónico",
        "iupac_name": "Nombre IUPAC (PubChem)",
        "stereogenic_elements": "Elementos Estereogénicos Detectados",
        "export_short_answer": "Generar Pregunta de Nomenclatura",
        "export_dragdrop": "Generar Pregunta de Estereoquímica (Drag & Drop)",
        "input_name_prompt": "Nombre IUPAC o referencia:",
        "success_generated": "¡Pregunta generada y añadida con éxito a la sesión!",
        "error_invalid_mol": "La estructura dibujada no es válida o está vacía.",
        "section_analysis": "Análisis Estructural",
        "section_export": "Exportación a Moodle"
    },
    "en": {
        "title": "Interactive Molecular Editor (Canvas) ✍️",
        "desc": "Draw chemical structures on the canvas, click 'Apply' inside the editor, and instantly get its analysis, IUPAC name, and Moodle export.",
        "canvas_label": "Chemical Drawing Canvas",
        "no_ketcher": "The `streamlit-ketcher` library is not installed.",
        "no_rdkit": "RDKit is not available in this environment.",
        "empty_canvas": "Draw a molecule on the canvas above and click 'Apply' to see its analysis.",
        "molecular_formula": "Molecular Formula",
        "molecular_weight": "Molecular Weight",
        "canonical_smiles": "Canonical SMILES",
        "iupac_name": "IUPAC Name (PubChem)",
        "stereogenic_elements": "Detected Stereogenic Elements",
        "export_short_answer": "Generate Nomenclature Question",
        "export_dragdrop": "Generate Stereochemistry Question (Drag & Drop)",
        "input_name_prompt": "IUPAC name or reference:",
        "success_generated": "Question successfully generated and added to session!",
        "error_invalid_mol": "The drawn structure is invalid or empty.",
        "section_analysis": "Structural Analysis",
        "section_export": "Moodle Export"
    }
}

def get_pubchem_iupac(smiles):
    """Consulta la API de PubChem para obtener el nombre IUPAC a partir del SMILES."""
    try:
        encoded_smiles = urllib.parse.quote(smiles)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded_smiles}/property/IUPACName/json"
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=3) as response:
            data = json.loads(response.read().decode())
            properties = data.get("PropertyTable", {}).get("Properties", [])
            if properties:
                return properties[0].get("IUPACName", None)
    except Exception:
        pass
    return None

def render_nueva_herramienta(texts_main):
    lang = st.session_state.get("lang", "en")
    t = TEXTS_EDITOR[lang]
    
    st.title(t["title"])
    st.write(t["desc"])
    
    if not RDKIT_AVAILABLE:
        st.warning(t["no_rdkit"])
        return
        
    if not KETCHER_AVAILABLE:
        st.warning(t["no_ketcher"])
        return

    if "nom_questions" not in st.session_state:
        st.session_state.nom_questions = []
    if "dd_questions" not in st.session_state:
        st.session_state.dd_questions = []

    st.markdown(f"### {t['canvas_label']}")
    
    # Render Ketcher Canvas
    smiles_drawn = st_ketcher("", height=480)
    
    if smiles_drawn and smiles_drawn.strip():
        mol = Chem.MolFromSmiles(smiles_drawn)
        if not mol:
            st.warning(t["error_invalid_mol"])
            return
            
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        formula = CalcMolFormula(mol)
        mw = Descriptors.MolWt(mol)
        
        # Search PubChem con caché
        cache_key = f"pubchem_name_{canonical_smiles}"
        if cache_key not in st.session_state:
            with st.spinner("Consultando nombre IUPAC en PubChem..." if lang=="es" else "Fetching IUPAC name from PubChem..."):
                fetched_name = get_pubchem_iupac(canonical_smiles)
                st.session_state[cache_key] = fetched_name if fetched_name else formula
        
        suggested_name = st.session_state[cache_key]
        
        # Refresh text box in session_state if molecule changes
        if "last_processed_smiles" not in st.session_state or st.session_state.last_processed_smiles != canonical_smiles:
            st.session_state.last_processed_smiles = canonical_smiles
            st.session_state.input_q_name_canvas = suggested_name
        
        chiral_centers = FindMolChiralCenters(mol, includeUnassigned=True)
        double_bonds = []
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                stereo = bond.GetStereo()
                if stereo in [Chem.rdchem.BondStereo.STEREOZ, Chem.rdchem.BondStereo.STEREOE]:
                    double_bonds.append(bond.GetIdx())

        st.markdown("---")
        
        # Create two paralel columns for analysis and export
        col_analysis, col_export = st.columns(2, gap="large")
        
        with col_analysis:
            st.markdown(f"#### 🔍 {t['section_analysis']}")
            st.markdown(f"**{t['canonical_smiles']}:** `{canonical_smiles}`")
            st.markdown(f"**{t['molecular_formula']}:** `{formula}`")
            st.markdown(f"**{t['molecular_weight']}:** `{mw:.2f} g/mol`")
            st.markdown(f"**{t['iupac_name']}:** *{suggested_name}*")
            
            stereo_desc = []
            for idx, cfg in chiral_centers:
                if cfg in ["R", "S"]:
                    stereo_desc.append(f"Atom #{idx+1}: {cfg}")
            for b_idx in double_bonds:
                b = mol.GetBondWithIdx(b_idx)
                stereo_desc.append(f"Bond #{b_idx} ({b.GetBeginAtomIdx()+1}={b.GetEndAtomIdx()+1}): {b.GetStereo().name}")
                
            if stereo_desc:
                st.markdown(f"**{t['stereogenic_elements']}:** " + ", ".join(stereo_desc))
            else:
                st.markdown(f"**{t['stereogenic_elements']}:** None / Achiral")
                
        with col_export:
            st.markdown(f"#### 📤 {t['section_export']}")
            
            # Use st.text_input directly connected to session_state to receive the updated IUPAC name
            q_name_input = st.text_input(t["input_name_prompt"], key="input_q_name_canvas")
                
            col_btn1, col_btn2 = st.columns(2)
            with col_btn1:
                submit_nom = st.button(
                    label=t["export_short_answer"], 
                    type="primary", 
                    icon=":material/edit_note:", 
                    use_container_width=True,
                    key="btn_export_nom_canvas"
                )
            with col_btn2:
                submit_dd = st.button(
                    label=t["export_dragdrop"], 
                    type="secondary", 
                    icon=":material/edit_note:", 
                    use_container_width=True,
                    key="btn_export_dd_canvas"
                )
            
            if submit_nom:
                opts = Chem.Draw.MolDrawOptions()
                opts.addAtomIndices = False
                img = Chem.Draw.MolToImage(mol, size=(400, 350), drawOptions=opts)
                buffered = io.BytesIO()
                img.save(buffered, format="PNG")
                img_b64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
                
                st.session_state.nom_questions.append({
                    "ref_name": q_name_input,
                    "valid_names": [q_name_input],
                    "smiles": canonical_smiles,
                    "image_b64": img_b64
                })
                st.success(t["success_generated"] + " (Nomenclature)")
                
            if submit_dd:
                if not chiral_centers and not double_bonds:
                    st.warning("La molécula no tiene elementos estereogénicos para Drag & Drop." if lang=="es" else "Molecule has no stereogenic elements for Drag & Drop.")
                else:
                    img_size = (500, 400)
                    drawer = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
                    opts = drawer.drawOptions()
                    opts.addAtomIndices = False
                    
                    highlight_atoms = [idx for idx, cfg in chiral_centers if cfg in ["R", "S"]]
                    drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, highlightBonds=double_bonds)
                    drawer.FinishDrawing()
                    
                    png_bytes = drawer.GetDrawingText()
                    img_b64 = base64.b64encode(png_bytes).decode("utf-8")
                    
                    elements = []
                    for idx, cfg in chiral_centers:
                        if cfg in ["R", "S"]:
                            pt = drawer.GetDrawCoords(idx)
                            elements.append({"index": idx+1, "label": cfg, "cx": float(pt.x), "cy": float(pt.y), "type": "atom"})
                            
                    for b_idx in double_bonds:
                        b = mol.GetBondWithIdx(b_idx)
                        stereo = b.GetStereo()
                        label_ze = "Z" if stereo == Chem.rdchem.BondStereo.STEREOZ else "E"
                        begin_idx = b.GetBeginAtomIdx()
                        end_idx = b.GetEndAtomIdx()
                        pt_b = drawer.GetDrawCoords(begin_idx)
                        pt_e = drawer.GetDrawCoords(end_idx)
                        mid_x = (float(pt_b.x) + float(pt_e.x)) / 2.0
                        mid_y = (float(pt_b.y) + float(pt_e.y)) / 2.0
                        elements.append({"index": begin_idx+1, "label": label_ze, "cx": mid_x, "cy": mid_y, "type": "bond"})
                        
                    st.session_state.dd_questions.append({
                        "name": q_name_input,
                        "image_b64": img_b64,
                        "elements": elements
                    })
                    st.success(t["success_generated"] + " (Drag & Drop)")
            else:
                    st.info(t["empty_canvas"])