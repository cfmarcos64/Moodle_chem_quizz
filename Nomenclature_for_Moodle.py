# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 09:05:18 2026

@author: Carlos Fernandez
"""

"""
Systematic Nomenclature Module for the Moodle Chemistry Suite.
Allows individual and bulk entry (Excel/CSV), clean molecular rendering,
alternative names support, field clearing on submit, Moodle XML export,
and detailed bulk upload error reporting.
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
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Bilingual interface texts
TEXTS_NOM = {
    "es": {
        "title": "Nomenclatura Sistemática de Moléculas 🏷️",
        "desc": "Genera preguntas Moodle de respuesta corta donde los alumnos deben identificar la nomenclatura sistemática IUPAC a partir de su estructura.",
        "tab_individual": "Entrada Individual",
        "tab_bulk": "Carga Masiva",
        "section_export": "Listado y Exportación de Preguntas",
        "compound_name": "Nombre base o referencia (ej: acetic acid)",
        "systematic_name": "Nombre sistemático IUPAC principal",
        "alt_name": "Nombre alternativo / sinónimo (opcional)",
        "add_question": "Validar y Añadir a la Lista",
        "bulk_uploader": "Sube un archivo Excel (.xlsx) o CSV",
        "bulk_note": "Columnas requeridas: **'compound'** y **'name_es'** (o 'name_en'). Opcional: **'alt_name'**.",
        "start_bulk": "Procesar Carga Masiva",
        "question_title": "Nomenclatura IUPAC: {}",
        "question_text": "¿Cuál es la nomenclatura sistemática IUPAC de esta molécula?",
        "export_xml": "Descargar XML Moodle (Respuesta Corta)",
        "clear_all": "Borrar Toda la Lista",
        "error_rdkit": "RDKit no está disponible en este entorno.",
        "error_smiles": "No se pudo encontrar un SMILES válido para '{}' en NCI CIR.",
        "error_required": "Por favor, completa los campos obligatorios.",
        "error_col_compound": "El archivo debe contener la columna 'compound'.",
        "success_added": "Pregunta añadida correctamente.",
        "bulk_success": "Carga masiva completada: **{}** éxitos, **{}** fallos.",
        "bulk_errors_title": "Ver detalle de errores en la carga masiva",
        "empty_list": "Aún no hay preguntas añadidas. Completa el formulario de la izquierda.",
        "last_added": "Última molécula añadida",
        "accepted_names": "Nombres aceptados"
    },
    "en": {
        "title": "Systematic Nomenclature of Molecules 🏷️",
        "desc": "Generate Moodle short-answer questions where students must identify the IUPAC systematic nomenclature from its structure.",
        "tab_individual": "Individual Entry",
        "tab_bulk": "Bulk Upload",
        "section_export": "Questions List & Export",
        "compound_name": "Base compound name or reference (e.g., acetic acid)",
        "systematic_name": "Main IUPAC systematic name",
        "alt_name": "Alternative name / synonym (optional)",
        "add_question": "Validate and Add to List",
        "bulk_uploader": "Upload Excel (.xlsx) or CSV file",
        "bulk_note": "Required columns: **'compound'** and **'name_en'** (or 'name_es'). Optional: **'alt_name'**.",
        "start_bulk": "Process Bulk Upload",
        "question_title": "IUPAC Nomenclature: {}",
        "question_text": "What is the IUPAC systematic nomenclature of this molecule?",
        "export_xml": "Download Moodle XML (Short Answer)",
        "clear_all": "Clear All List",
        "error_rdkit": "RDKit is not available in this environment.",
        "error_smiles": "Could not find a valid SMILES for '{}' via NCI CIR.",
        "error_required": "Please fill in required fields.",
        "error_col_compound": "The file must contain a 'compound' column.",
        "success_added": "Question successfully added.",
        "bulk_success": "Bulk processing complete: **{}** successes, **{}** failures.",
        "bulk_errors_title": "View bulk upload error details",
        "empty_list": "No questions added yet. Fill out the form on the left.",
        "last_added": "Last added molecule",
        "accepted_names": "Accepted names"
    }
}

def resolver_smiles(input_str):
    """Retrieves canonical SMILES from a name (via NCI CIR) or validates direct SMILES."""
    input_limpio = str(input_str).strip()
    if not input_limpio:
        return None
    
    if RDKIT_AVAILABLE:
        mol_test = Chem.MolFromSmiles(input_limpio)
        if mol_test:
            return Chem.MolToSmiles(mol_test, canonical=True)
            
    try:
        encoded = requests.utils.quote(input_limpio)
        url = f"http://cactus.nci.nih.gov/chemical/structure/{encoded}/smiles"
        res = requests.get(url, timeout=10)
        if res.status_code == 200:
            smiles = res.text.strip().split('\n')[0]
            if smiles and "Error" not in smiles:
                if RDKIT_AVAILABLE:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        return Chem.MolToSmiles(mol, canonical=True)
                return smiles
    except Exception:
        pass
    return None

def generar_imagen_limpia_b64(smiles):
    """Generates a clean Base64 PNG image of the molecule (without atom indices)."""
    if not RDKIT_AVAILABLE:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            opts = Draw.MolDrawOptions()
            opts.addAtomIndices = False
            img = Draw.MolToImage(mol, size=(400, 350), drawOptions=opts)
            buffered = io.BytesIO()
            img.save(buffered, format="PNG")
            return base64.b64encode(buffered.getvalue()).decode("utf-8")
    except Exception:
        pass
    return None

def generar_short_answer_xml(questions_list, lang):
    """Generates a Moodle-compatible XML file for short-answer questions using standard HTML text formatting."""
    t = TEXTS_NOM[lang]
    root = ET.Element("quiz")
     
    for item in questions_list:
        ref_name = item["ref_name"]
        names_validos = item["valid_names"]
        b64_img = item["image_b64"]
        
        question = ET.Element("question", type="shortanswer")
        
        # Question name (Corregido para evitar la duplicación)
        q_name = ET.SubElement(question, "name")
        ET.SubElement(q_name, "text").text = t["question_title"].format(ref_name)
        
        # Question text with embedded image using standard HTML text nodes
        q_text = ET.SubElement(question, "questiontext", format="html")
        img_tag = f'<p><img src="data:image/png;base64,{b64_img}" alt="{ref_name}" /></p>' if b64_img else ""
        ET.SubElement(q_text, "text").text = f"<p>{t['question_text']}</p>{img_tag}"
        
        # Configure valid answers (all with 100% score)
        for n_val in names_validos:
            answer = ET.SubElement(question, "answer", fraction="100", format="moodle_auto_format")
            ET.SubElement(answer, "text").text = n_val
            feedback = ET.SubElement(answer, "feedback", format="html")
            ET.SubElement(feedback, "text").text = "¡Correcto!" if lang == "es" else "Correct!"
            
        root.append(question)
        
    buffer = io.BytesIO()
    buffer.write(ET.tostring(root, encoding="utf-8", xml_declaration=True))
    return buffer.getvalue()

def render_nueva_herramienta(texts_main):
    """Main function that renders the nomenclature module in Streamlit."""
    lang = st.session_state.get("lang", "en")
    t = TEXTS_NOM[lang]
    
    st.title(t["title"])
    st.write(t["desc"])
    
    if not RDKIT_AVAILABLE:
        st.warning(t["error_rdkit"])
        return

    # Initialize state for this module
    if "nom_questions" not in st.session_state:
        st.session_state.nom_questions = []
        
    if "last_nom_preview" not in st.session_state:
        st.session_state.last_nom_preview = None

    # --- TWO-COLUMN LAYOUT ---
    col_input, col_list = st.columns([1.1, 0.9], gap="large")

    # ================= LEFT COLUMN: INPUT (INDIVIDUAL & BULK) =================
    with col_input:
        tab_ind, tab_bulk = st.tabs([t["tab_individual"], t["tab_bulk"]])

        with tab_ind:
            with st.form("form_individual_nom", clear_on_submit=True):
                c_name = st.text_input(t["compound_name"])
                sys_name = st.text_input(t["systematic_name"])
                alt_name = st.text_input(t["alt_name"])
                    
                submitted = st.form_submit_button(t["add_question"], type="primary", icon=":material/add_task:")
                
                if submitted:
                    if not c_name.strip() or not sys_name.strip():
                        st.error(t["error_required"])
                    else:
                        with st.spinner("Searching structure and generating image..." if lang=="en" else "Buscando estructura y generando imagen..."):
                            smiles = resolver_smiles(c_name)
                            if not smiles:
                                st.error(t["error_smiles"].format(c_name))
                            else:
                                img_b64 = generar_imagen_limpia_b64(smiles)
                                validos = [sys_name.strip()]
                                if alt_name.strip():
                                    validos.append(alt_name.strip())
                                    
                                st.session_state.nom_questions.append({
                                    "ref_name": c_name.strip(),
                                    "valid_names": validos,
                                    "smiles": smiles,
                                    "image_b64": img_b64
                                })
                                
                                # Store preview to show outside the form
                                st.session_state.last_nom_preview = {
                                    "ref_name": c_name.strip(),
                                    "valid_names": validos,
                                    "image_b64": img_b64
                                }
                                st.success(t["success_added"])

            # Show preview of the last added compound outside the form
            if st.session_state.last_nom_preview:
                prev = st.session_state.last_nom_preview
                nombres_str = " / ".join([f"`{n}`" for n in prev['valid_names']])
                st.markdown(f"**{t['last_added']}:** `{prev['ref_name']}`<br>{t['accepted_names']}: {nombres_str}", unsafe_allow_html=True)
                if prev["image_b64"]:
                    st.image(f"data:image/png;base64,{prev['image_b64']}", width=250, caption=prev["ref_name"])

        with tab_bulk:
            st.info(t["bulk_note"])
            uploaded_file = st.file_uploader(t["bulk_uploader"], type=["xlsx", "csv"], key="bulk_nom_file")
            
            if uploaded_file and st.button(t["start_bulk"]):
                try:
                    if uploaded_file.name.endswith(".csv"):
                        df = pd.read_csv(uploaded_file)
                    else:
                        df = pd.read_excel(uploaded_file)
                        
                    if not "compound" in df.columns:
                        st.error(t["error_col_compound"])
                    else:
                        success_count = 0
                        fail_count = 0
                        errores_detalle = []
                        
                        progress_bar = st.progress(0, text="Processing batch..." if lang=="en" else "Procesando lote...")
                        total = len(df)
                        
                        for idx, row in df.iterrows():
                            comp_input = str(row["compound"]).strip() if pd.notna(row["compound"]) else ""
                            col_target = "name_es" if lang == "es" and "name_es" in df.columns else ("name_en" if "name_en" in df.columns else (df.columns[1] if len(df.columns) > 1 else None))
                            
                            n_val = str(row[col_target]).strip() if col_target and col_target in df.columns and pd.notna(row[col_target]) else ""
                            alt_val = str(row["alt_name"]).strip() if "alt_name" in df.columns and pd.notna(row["alt_name"]) else ""
                            
                            progress_bar.progress((idx + 1) / total, text=f"Processing: {comp_input}")
                            
                            if not comp_input:
                                fail_count += 1
                                errores_detalle.append({"row": idx + 1, "name": "N/A", "cause": "Missing compound name"})
                                continue
                                
                            if not n_val:
                                fail_count += 1
                                errores_detalle.append({"row": idx + 1, "name": comp_input, "cause": "Missing systematic name (name_en/name_es)"})
                                continue
                                
                            smiles = resolver_smiles(comp_input)
                            if not smiles:
                                fail_count += 1
                                errores_detalle.append({"row": idx + 1, "name": comp_input, "cause": "SMILES not found in NCI CIR"})
                                continue
                                
                            img_b64 = generar_imagen_limpia_b64(smiles)
                            validos = [n_val]
                            if alt_val:
                                validos.append(alt_val)
                                
                            st.session_state.nom_questions.append({
                                "ref_name": comp_input,
                                "valid_names": validos,
                                "smiles": smiles,
                                "image_b64": img_b64
                            })
                            success_count += 1
                                
                        progress_bar.empty()
                        st.success(t["bulk_success"].format(success_count, fail_count))
                        
                        # Show error details if any
                        if errores_detalle:
                            with st.expander(t["bulk_errors_title"]):
                                df_err = pd.DataFrame(errores_detalle)
                                st.dataframe(df_err, use_container_width=True)
                                
                except Exception as e:
                    st.error(f"Error processing file: {e}")

    # ================= RIGHT COLUMN: LIST & EXPORT =================
    with col_list:
        st.subheader(t["section_export"])
        
        if st.session_state.nom_questions:
            st.write(f"Total questions in list: **{len(st.session_state.nom_questions)}**" if lang=="en" else f"Total de preguntas en la lista: **{len(st.session_state.nom_questions)}**")
            
            # --- VERTICALLY ALIGNED MAIN BUTTONS ---
            xml_data = generar_short_answer_xml(st.session_state.nom_questions, lang)
            
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
                file_name="moodle_nomenclature.xml",
                icon=":material/download:",
                mime="application/xml",
                use_container_width=True
            )

            if st.button(t["clear_all"], type="secondary", icon=":material/delete_sweep:", use_container_width=True):
                st.session_state.nom_questions = []
                st.session_state.last_nom_preview = None
                st.rerun()
            
            st.divider()
                
            # Scrollable container for questions list
            container_list = st.container(height=420)
            with container_list:
                for i, q in enumerate(st.session_state.nom_questions):
                    cols_item = st.columns([1, 2.5, 0.6])
                    with cols_item[0]:
                        if q['image_b64']:
                            st.image(f"data:image/png;base64,{q['image_b64']}", width=100)
                    with cols_item[1]:
                        nombres_str = " / ".join([f"`{n}`" for n in q['valid_names']])
                        label_accepted = "Accepted" if lang == "en" else "Aceptados"
                        st.markdown(f"**#{i+1}** ({q['ref_name']})<br>{label_accepted}: {nombres_str}", unsafe_allow_html=True)
                    with cols_item[2]:
                        if st.button("🗑️", key=f"del_nom_{i}"):
                            st.session_state.nom_questions.pop(i)
                            st.rerun()
                    st.divider()
        else:
            st.info(t["empty_list"])
