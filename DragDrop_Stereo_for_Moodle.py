# -*- coding: utf-8 -*-
"""
Interactive Stereochemistry Module (Drag & Drop - Moodle ddimageortext)
Automatically extracts exact pixel coordinates of stereocenters and bond midpoints using RDKit MolDraw2DCairo.
"""

import streamlit as st
import xml.etree.ElementTree as ET
import pandas as pd
import requests
import io
import base64

# --- Dependency Verification (RDKit) ---
try:
    from rdkit import Chem
    from rdkit.Chem import FindMolChiralCenters
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Helper for CDATA in ElementTree
class CDATA(str):
    pass

def _cdata_serializer(current, parent, index, *args, **kwargs):
    if isinstance(current, CDATA):
        current.tag = 'text'
        element = ET.Element('text')
        element.text = current
        return element
    return None

# Bilingual interface texts
TEXTS_DD = {
    "es": {
        "title": "Estereoquímica Interactiva (Drag & Drop) 🎯",
        "desc": "Detecta elementos estereogénicos (R/S y Z/E), calcula sus coordenadas exactas en la imagen y genera preguntas de Moodle de arrastrar y soltar.",
        "tab_individual": "Entrada Individual",
        "tab_bulk": "Carga Masiva",
        "section_export": "Listado y Exportación de Preguntas",
        "compound_name": "Nombre de la molécula en inglés (ej: (S)-butan-2-ol, (E)-but-2-ene)",
        "add_question": "Analizar y Añadir Pregunta Interactiva",
        "bulk_uploader": "Sube un archivo Excel (.xlsx) o CSV",
        "bulk_note": "El archivo debe contener al menos una columna llamada **'name'** o **'nombre'**.",
        "start_bulk": "Procesar Carga Masiva",
        "export_xml": "Descargar XML Moodle (Drag & Drop)",
        "clear_all": "Borrar Toda la Lista",
        "error_rdkit": "RDKit no está disponible en este entorno.",
        "error_smiles": "No se pudo encontrar un SMILES válido para '{}' en NCI CIR.",
        "error_chiral": "La molécula no presenta elementos estereogénicos detectables.",
        "success_added": "Pregunta interactiva añadida con éxito.",
        "empty_list": "Aún no hay preguntas añadidas.",
        "last_added": "Última molécula añadida",
        "elements_detected": "Elementos estereogénicos detectados"
    },
    "en": {
        "title": "Interactive Stereochemistry (Drag & Drop) 🎯",
        "desc": "Detects stereogenic elements (R/S & Z/E), computes their exact image coordinates, and generates Moodle drag-and-drop questions.",
        "tab_individual": "Individual Entry",
        "tab_bulk": "Bulk Upload",
        "section_export": "Questions List & Export",
        "compound_name": "Molecule name in English (e.g., (S)-butan-2-ol, (E)-but-2-ene)",
        "add_question": "Analyze and Add Interactive Question",
        "bulk_uploader": "Upload Excel (.xlsx) or CSV file",
        "bulk_note": "The file must contain a column named **'name'** or **'nombre'**.",
        "start_bulk": "Process Bulk Upload",
        "export_xml": "Download Moodle XML (Drag & Drop)",
        "clear_all": "Clear All Questions",
        "error_rdkit": "RDKit is not available in this environment.",
        "error_smiles": "Could not find a valid SMILES for '{}' via NCI CIR.",
        "error_chiral": "The molecule has no detectable stereogenic elements.",
        "success_added": "Interactive question successfully added.",
        "empty_list": "No questions added yet.",
        "last_added": "Last added molecule",
        "elements_detected": "Detected stereogenic elements"
    }
}

def buscar_smiles_nci(nombre_molecula):
    """Searches for the canonical SMILES of a molecule using the NCI CIR API."""
    nombre_limpio = nombre_molecula.strip()
    url = f"http://cactus.nci.nih.gov/chemical/structure/{requests.utils.quote(nombre_limpio)}/smiles"
    try:
        respuesta = requests.get(url, timeout=8)
        if respuesta.status_code == 200:
            smiles = respuesta.text.strip().split('\n')[0]
            if smiles and "Error" not in smiles:
                return smiles
    except Exception:
        pass
    return None

def analizar_drag_drop_mol(nombre_molecula):
    """
    Analyzes stereochemistry and extracts exact pixel coordinates of stereocenters
    and bond midpoints using RDKit MolDraw2DCairo.
    """
    if not RDKIT_AVAILABLE:
        return None, [], "RDKit not available"

    smiles = buscar_smiles_nci(nombre_molecula)
    if not smiles:
        return None, [], f"SMILES not found for {nombre_molecula}"

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None, [], f"Error parsing SMILES: {smiles}"

    try:
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    except Exception as e:
        return None, [], f"Error assigning stereochemistry: {e}"

    img_size = (500, 400)
    
    elementos_estereo = []
    indices_resaltar_atomos = []
    indices_resaltar_enlaces = []

    # 1. Chiral centers (R/S)
    centros_quirales = FindMolChiralCenters(mol, includeUnassigned=True)
    for idx, config in centros_quirales:
        if config in ["R", "S"]:
            indices_resaltar_atomos.append(idx)

    # 2. Geometric double bonds (Z/E)
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            stereo = bond.GetStereo()
            if stereo in [Chem.rdchem.BondStereo.STEREOZ, Chem.rdchem.BondStereo.STEREOE]:
                indices_resaltar_enlaces.append(bond.GetIdx())

    if not centros_quirales and not indices_resaltar_enlaces:
        return None, [], "No detectable stereogenic elements"

    # Use MolDraw2DCairo to draw image and get exact coordinates from the same canvas
    drawer = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
    opts = drawer.drawOptions()
    opts.addAtomIndices = False
    
    drawer.DrawMolecule(
        mol,
        highlightAtoms=indices_resaltar_atomos,
        highlightBonds=indices_resaltar_enlaces
    )
    drawer.FinishDrawing()
    
    png_bytes = drawer.GetDrawingText()
    img_b64 = base64.b64encode(png_bytes).decode("utf-8")

    # Extract exact coordinates for chiral centers (R/S)
    for idx, config in centros_quirales:
        if config in ["R", "S"]:
            iupac_idx = idx + 1
            pt = drawer.GetDrawCoords(idx)
            elementos_estereo.append({
                "index": iupac_idx,
                "label": config,
                "cx": float(pt.x),
                "cy": float(pt.y),
                "type": "atom"
            })

    # Extract midpoint coordinates for double bonds (Z/E) between the two atoms forming the bond
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            stereo = bond.GetStereo()
            if stereo in [Chem.rdchem.BondStereo.STEREOZ, Chem.rdchem.BondStereo.STEREOE]:
                label_ze = "Z" if stereo == Chem.rdchem.BondStereo.STEREOZ else "E"
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                iupac_idx = begin_idx + 1
                
                pt_begin = drawer.GetDrawCoords(begin_idx)
                pt_end = drawer.GetDrawCoords(end_idx)
                
                mid_x = (float(pt_begin.x) + float(pt_end.x)) / 2.0
                mid_y = (float(pt_begin.y) + float(pt_end.y)) / 2.0
                
                elementos_estereo.append({
                    "index": iupac_idx,
                    "label": label_ze,
                    "cx": mid_x,
                    "cy": mid_y,
                    "type": "bond"
                })

    if not elementos_estereo:
        return None, [], "No detectable stereogenic elements"

    return img_b64, elementos_estereo, None

def generar_xml_ddimageortext(questions_list, lang):
    """Generates Moodle XML for ddimageortext matching database requirements."""
    quiz = ET.Element("quiz")
    
    for item in questions_list:
        name = item["name"]
        b64_img = item["image_b64"]
        elements = item["elements"]
        
        question = ET.Element("question", type="ddimageortext")
        
        # Question name
        q_name = ET.SubElement(question, "name")
        q_text_name = ET.SubElement(q_name, "text")
        q_text_name.text = f"Drag & Drop Stereo: {name}"
        
        # Question text with embedded image using pluginfile syntax
        q_text = ET.SubElement(question, "questiontext", format="html")
        q_text_content = ET.SubElement(q_text, "text")
        
        prompt = "Arrastra las etiquetas estereoquímicas correctas sobre los elementos señalados de la molécula:" if lang == "es" else "Drag the correct stereochemical labels onto the specified elements of the molecule:"
        img_tag = f'<p><img src="@@PLUGINFILE@@/molecule.png" alt="{name}" /></p>'
        q_text_content.text = f"<p>{prompt}</p>{img_tag}"
        
        # General feedback and grading parameters required by ddimageortext
        feedback = ET.SubElement(question, "generalfeedback", format="html")
        ET.SubElement(feedback, "text").text = ""
        
        ET.SubElement(question, "defaultgrade").text = "1.0000000"
        ET.SubElement(question, "penalty").text = "0.3333333"
        ET.SubElement(question, "hidden").text = "0"
        ET.SubElement(question, "idnumber").text = ""
        ET.SubElement(question, "correctfeedback", format="html").append(ET.Element("text"))
        ET.SubElement(question, "incorrectfeedback", format="html").append(ET.Element("text"))
        ET.SubElement(question, "shownumcorrect").text = ""

        # File node for background image in Moodle ddimageortext
        file_node = ET.SubElement(question, "file", name="molecule.png", encoding="base64")
        file_node.text = b64_img
        
        # Draggable items: Always include R and S if there are chiral centers, Z and E if there are double bonds.
        has_chiral = any(el["label"] in ["R", "S"] for el in elements)
        has_geometric = any(el["label"] in ["Z", "E"] for el in elements)
        
        drags_set = []
        if has_chiral:
            drags_set.extend(["R", "S"])
        if has_geometric:
            drags_set.extend(["Z", "E"])
        if not drags_set:
            drags_set = sorted(list({el["label"] for el in elements}))

        drag_mapping = {}
        for i, drag_text in enumerate(drags_set):
            drag_no = i + 1
            drag_mapping[drag_text] = drag_no
            drag = ET.SubElement(question, "drag")
            ET.SubElement(drag, "no").text = str(drag_no)
            ET.SubElement(drag, "group").text = "1"
            ET.SubElement(drag, "infinite").text = "1"
            drag_text_node = ET.SubElement(drag, "text")
            drag_text_node.text = drag_text
            
        # Drop zones with full compatibility (xleft, ytop, shape, coords)
        for i, el in enumerate(elements):
            drop = ET.SubElement(question, "drop")
            ET.SubElement(drop, "no").text = str(i + 1)
            ET.SubElement(drop, "choice").text = str(drag_mapping[el["label"]])
            
            cx = el["cx"]
            cy = el["cy"]
            left = max(0, int(cx - 20))
            top = max(0, int(cy - 20))
            width = 40
            height = 40
            
            ET.SubElement(drop, "xleft").text = str(left)
            ET.SubElement(drop, "ytop").text = str(top)
            ET.SubElement(drop, "shape").text = "rectangle"
            ET.SubElement(drop, "coords").text = f"{left},{top},{width},{height}"
            
        quiz.append(question)
        
    tree = ET.ElementTree(quiz)
    if hasattr(ET, "indent"):
        ET.indent(tree, space="    ")
        
    buffer = io.BytesIO()
    tree.write(buffer, encoding="utf-8", xml_declaration=True)
    return buffer.getvalue()

def render_nueva_herramienta(texts_main):
    """Main function rendering the Drag & Drop module in Streamlit."""
    lang = st.session_state.get("lang", "en")
    t = TEXTS_DD[lang]
    
    st.title(t["title"])
    st.write(t["desc"])
    
    if not RDKIT_AVAILABLE:
        st.warning(t["error_rdkit"])
        return

    if "dd_questions" not in st.session_state:
        st.session_state.dd_questions = []
        
    if "last_dd_preview" not in st.session_state:
        st.session_state.last_dd_preview = None

    col_input, col_list = st.columns([1.1, 0.9], gap="large")

    with col_input:
        tab_ind, tab_bulk = st.tabs([t["tab_individual"], t["tab_bulk"]])

        with tab_ind:
            with st.form("form_individual_dd", clear_on_submit=True):
                c_name = st.text_input(t["compound_name"])
                submitted = st.form_submit_button(t["add_question"], type="primary", icon=":material/extension:")
                
                if submitted:
                    if not c_name.strip():
                        st.error("Por favor, introduce el nombre de la molécula." if lang=="es" else "Please enter a molecule name.")
                    else:
                        with st.spinner("Analizando geometría y calculando zonas de caída..." if lang=="es" else "Analyzing geometry and calculating drop zones..."):
                            img_b64, elements, err = analizar_drag_drop_mol(c_name)
                            if err:
                                st.warning(f"Aviso: {err}" if lang=="es" else f"Notice: {err}")
                            elif not elements:
                                st.warning(t["error_chiral"])
                            else:
                                st.session_state.dd_questions.append({
                                    "name": c_name.strip(),
                                    "image_b64": img_b64,
                                    "elements": elements
                                })
                                st.session_state.last_dd_preview = {
                                    "name": c_name.strip(),
                                    "image_b64": img_b64,
                                    "elements": elements
                                }
                                st.success(t["success_added"])

            if st.session_state.last_dd_preview:
                prev = st.session_state.last_dd_preview
                st.markdown(f"**{t['last_added']}:** `{prev['name']}`")
                elems_str = ", ".join([f"Pos {e['index']}: {e['label']}" for e in prev['elements']])
                st.markdown(f"**{t['elements_detected']}:** {elems_str}")
                if prev["image_b64"]:
                    st.image(f"data:image/png;base64,{prev['image_b64']}", width=250, caption=prev["name"])

        with tab_bulk:
            st.info(t["bulk_note"])
            uploaded_file = st.file_uploader(t["bulk_uploader"], type=["xlsx", "csv"], key="uploader_dd")
            
            if uploaded_file and st.button(t["start_bulk"]):
                try:
                    if uploaded_file.name.endswith(".csv"):
                        df = pd.read_csv(uploaded_file)
                    else:
                        df = pd.read_excel(uploaded_file)
                        
                    col_name = 'name' if 'name' in df.columns else ('nombre' if 'nombre' in df.columns else None)
                    if not col_name:
                        st.error("El archivo debe contener una columna 'name'." if lang=="es" else "File must contain a 'name' column.")
                    else:
                        total_rows = len(df)
                        success = 0
                        error_details = []
                        
                        progress_bar = st.progress(0)
                        status_text = st.empty()
                        
                        for idx, row in df.iterrows():
                            row_num = idx + 2 # Considerando cabecera en fila 1 (Excel index 0 -> row 2)
                            c_name = str(row[col_name]).strip()
                            
                            # Actualizar barra de progreso y texto de estado
                            pct = float(idx + 1) / total_rows
                            progress_bar.progress(pct)
                            status_text.text(f"Procesando fila {row_num}/{total_rows}: `{c_name}`..." if lang=="es" else f"Processing row {row_num}/{total_rows}: `{c_name}`...")
                            
                            img_b64, elements, err = analizar_drag_drop_mol(c_name)
                            if err or not elements:
                                reason = err if err else ("Sin elementos estereogénicos" if lang=="es" else "No stereogenic elements")
                                error_details.append(f"Fila {row_num} (`{c_name}`): {reason}")
                            else:
                                st.session_state.dd_questions.append({
                                    "name": c_name,
                                    "image_b64": img_b64,
                                    "elements": elements
                                })
                                success += 1
                        
                        status_text.empty()
                        progress_bar.empty()
                        
                        st.success(f"Carga masiva completada: **{success}** añadidas con éxito, **{len(error_details)}** omitidas." if lang=="es" else f"Bulk complete: **{success}** successfully added, **{len(error_details)}** skipped.")
                        
                        if error_details:
                            with st.expander(f"Ver detalles de líneas con errores u omitidas ({len(error_details)})" if lang=="es" else f"View skipped/error rows detail ({len(error_details)})"):
                                for detail in error_details:
                                    st.markdown(f"- {detail}")
                                    
                except Exception as e:
                    st.error(f"Error procesando archivo: {e}")

    with col_list:
        st.subheader(t["section_export"])
        if st.session_state.dd_questions:
            st.write(f"Total: **{len(st.session_state.dd_questions)}**")
            
            xml_data = generar_xml_ddimageortext(st.session_state.dd_questions, lang)
            
            st.markdown("""
                <style>
                div.stDownloadButton > button:first-child {
                    background-color: #28a745 !important;
                    color: white !important;
                    font-weight: bold !important;
                    border: none !important;
                    width: 100% !important;
                }
                div.stDownloadButton > button:first-child:hover {
                    background-color: #218838 !important;
                }
                </style>
            """, unsafe_allow_html=True)
            
            st.download_button(
                label=t["export_xml"],
                data=xml_data,
                file_name="moodle_dragdrop_stereochemistry.xml",
                icon=":material/download:",
                mime="application/xml",
                use_container_width=True
            )

            if st.button(t["clear_all"], type="secondary", icon=":material/delete_sweep:", use_container_width=True):
                st.session_state.dd_questions = []
                st.session_state.last_dd_preview = None
                st.rerun()
                
            st.divider()
            
            container_list = st.container(height=420)
            with container_list:
                for i, q in enumerate(st.session_state.dd_questions):
                    cols = st.columns([1, 2.5, 0.6])
                    with cols[0]:
                        if q['image_b64']:
                            st.image(f"data:image/png;base64,{q['image_b64']}", width=100)
                    with cols[1]:
                        elems = ", ".join([f"#{e['index']}:{e['label']}" for e in q['elements']])
                        st.markdown(f"**#{i+1}** (`{q['name']}`)<br>Zonas: {elems}", unsafe_allow_html=True)
                    with cols[2]:
                        if st.button("🗑️", key=f"del_dd_{i}"):
                            st.session_state.dd_questions.pop(i)
                            st.rerun()
                    st.divider()
        else:
            st.info(t["empty_list"])
