# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:51:39 2026

@author: Carlos
"""

"""
Módulo de Determinación de Grupo Funcional para la Suite de Moodle.
Permite entrada individual y masiva (Excel), renderizado de moléculas con RDKit,
pestañas separadas, diseño en dos columnas, limpieza automática de campos e imágenes en el listado.
"""

import streamlit as st
import xml.etree.ElementTree as ET
import pandas as pd
import requests
import io
import base64

# --- Verificación de Dependencias (RDKit) ---
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Textos bilingües para la interfaz
TEXTS_DG = {
    "es": {
        "title": "Determinación de Grupo Funcional 🧪",
        "desc": "Genera preguntas Moodle de respuesta corta sobre el grupo funcional principal de una molécula.",
        "tab_individual": "Entrada Individual",
        "tab_bulk": "Carga Masiva",
        "section_export": "Listado y Exportación de Preguntas",
        "compound_name": "Nombre de la molécula (en inglés, ej: acetic acid)",
        "functional_group": "Nombre del grupo funcional principal (respuesta esperada)",
        "add_question": "Validar y Añadir a la Lista",
        "bulk_uploader": "Sube un archivo Excel (.xlsx) o CSV",
        "bulk_note": "El archivo debe contener las columnas: **'name'** y **'functional_group'**.",
        "start_bulk": "Procesar Carga Masiva",
        "question_title": "Identificación de Grupo Funcional: {}",
        "question_text": "¿Cuál es el grupo funcional principal de esta molécula?",
        "export_xml": "Descargar XML Moodle (Respuesta Corta)",
        "clear_all": "Borrar Toda la Lista",
        "error_rdkit": "RDKit no está disponible en este entorno para generar las imágenes moleculares.",
        "error_smiles": "No se pudo encontrar un SMILES válido para el compuesto '{}' en NCI CIR.",
        "success_added": "Pregunta añadida correctamente.",
        "bulk_success": "Carga masiva completada: **{}** éxitos, **{}** fallos.",
        "empty_list": "Aún no hay preguntas añadidas. Completa el formulario de la izquierda."
    },
    "en": {
        "title": "Functional Group Determination 🧪",
        "desc": "Generate Moodle short-answer questions regarding the main functional group of a molecule.",
        "tab_individual": "Individual Entry",
        "tab_bulk": "Bulk Upload",
        "section_export": "Questions List & Export",
        "compound_name": "Molecule Name (in English, e.g., acetic acid)",
        "functional_group": "Main Functional Group (expected answer)",
        "add_question": "Validate and Add to List",
        "bulk_uploader": "Upload Excel (.xlsx) or CSV file",
        "bulk_note": "The file must contain columns: **'name'** and **'functional_group'**.",
        "start_bulk": "Process Bulk Upload",
        "question_title": "Functional Group Identification: {}",
        "question_text": "What is the main functional group of this molecule?",
        "export_xml": "Download Moodle XML (Short Answer)",
        "clear_all": "Clear All List",
        "error_rdkit": "RDKit is not available in this environment to generate molecular images.",
        "error_smiles": "Could not find a valid SMILES for compound '{}' via NCI CIR.",
        "success_added": "Question successfully added.",
        "bulk_success": "Bulk processing complete: **{}** successes, **{}** failures.",
        "empty_list": "No questions added yet. Fill out the form on the left."
    }
}

def name_to_smiles_simple(compound_name):
    """Busca el SMILES de una molécula usando la API de NCI CIR."""
    try:
        encoded = requests.utils.quote(compound_name.strip())
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

def smiles_to_png_base64(smiles):
    """Genera una imagen PNG en Base64 de la molécula usando RDKit de forma limpia."""
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

def generate_short_answer_xml(questions_list, lang):
    """Genera el XML para Moodle de tipo 'shortanswer' incluyendo imágenes integradas en Base64."""
    quiz = ET.Element("quiz")
    
    for item in questions_list:
        name = item["name"]
        fg = item["functional_group"]
        b64_img = item["image_b64"]
        
        question = ET.Element("question", type="shortanswer")
        
        # Nombre de la pregunta
        q_name = ET.SubElement(question, "name")
        q_text_name = ET.SubElement(q_name, "text")
        q_text_name.text = TEXTS_DG[lang]["question_title"].format(name)
        
        # Texto del enunciado con imagen embebida en base64
        q_text = ET.SubElement(question, "questiontext", format="html")
        q_text_content = ET.SubElement(q_text, "text")
        
        base_prompt = TEXTS_DG[lang]["question_text"]
        img_tag = f'<p><img src="data:image/png;base64,{b64_img}" alt="{name}" /></p>' if b64_img else ""
        q_text_content.text = f"<![CDATA[<p>{base_prompt}</p>{img_tag}]]>"
        
        # Opciones de respuesta corta
        answer = ET.SubElement(question, "answer", fraction="100", format="moodle_auto_format")
        ans_text = ET.SubElement(answer, "text")
        ans_text.text = fg
        
        feedback = ET.SubElement(answer, "feedback", format="html")
        ET.SubElement(feedback, "text").text = "¡Correcto!" if lang == "es" else "Correct!"
        
        quiz.append(question)
        
    tree = ET.ElementTree(quiz)
    ET.indent(tree, space="  ")
    return ET.tostring(tree.getroot(), encoding="utf-8", xml_declaration=True)

def render_nueva_herramienta(texts_main):
    """Función principal que renderiza el módulo en Streamlit."""
    lang = st.session_state.get("lang", "en")
    t = TEXTS_DG[lang]
    
    st.title(t["title"])
    st.write(t["desc"])
    
    if not RDKIT_AVAILABLE:
        st.warning(t["error_rdkit"])
        return

    # Inicializar estado para este módulo si no existe
    if "dg_questions" not in st.session_state:
        st.session_state.dg_questions = []
        
    if "last_fg_preview" not in st.session_state:
        st.session_state.last_fg_preview = None

    # --- DISEÑO DE DOS COLUMNAS ---
    col_input, col_list = st.columns([1.1, 0.9], gap="large")

    # ================= COLUMNA IZQUIERDA: ENTRADA (INDIVIDUAL Y MASIVA) =================
    with col_input:
        tab_ind, tab_bulk = st.tabs([t["tab_individual"], t["tab_bulk"]])

        with tab_ind:
            with st.form("form_individual_dg", clear_on_submit=True):
                c_name = st.text_input(t["compound_name"])
                f_group = st.text_input(t["functional_group"])
                submitted = st.form_submit_button(t["add_question"], type="primary", icon=":material/add_task:")
                
                if submitted:
                    if not c_name.strip() or not f_group.strip():
                        st.error("Por favor, completa ambos campos." if lang == "es" else "Please fill in all fields.")
                    else:
                        with st.spinner("Buscando estructura y generando imagen..."):
                            smiles = name_to_smiles_simple(c_name)
                            if not smiles:
                                st.error(t["error_smiles"].format(c_name))
                            else:
                                img_b64 = smiles_to_png_base64(smiles)
                                st.session_state.dg_questions.append({
                                    "name": c_name.strip(),
                                    "functional_group": f_group.strip(),
                                    "smiles": smiles,
                                    "image_b64": img_b64
                                })
                                # Guardar la vista previa para mostrarla fuera del formulario
                                st.session_state.last_fg_preview = {
                                    "name": c_name.strip(),
                                    "functional_group": f_group.strip(),
                                    "image_b64": img_b64
                                }
                                st.success(t["success_added"])

            # Mostrar la previsualización del último compuesto añadido fuera del formulario
            if st.session_state.last_fg_preview:
                prev = st.session_state.last_fg_preview
                st.markdown(f"**Última molécula añadida:** `{prev['name']}`<br>Grupo: **{prev['functional_group']}**", unsafe_allow_html=True)
                if prev["image_b64"]:
                    st.image(f"data:image/png;base64,{prev['image_b64']}", width=250, caption=prev["name"])

        with tab_bulk:
            st.info(t["bulk_note"])
            uploaded_file = st.file_uploader(t["bulk_uploader"], type=["xlsx", "csv"], key="bulk_dg_file")
            
            if uploaded_file and st.button(t["start_bulk"]):
                try:
                    if uploaded_file.name.endswith(".csv"):
                        df = pd.read_csv(uploaded_file)
                    else:
                        df = pd.read_excel(uploaded_file)
                        
                    if not {"name", "functional_group"}.issubset(df.columns):
                        st.error("El archivo debe contener las columnas exactas: 'name' y 'functional_group'.")
                    else:
                        success_count = 0
                        fail_count = 0
                        progress_bar = st.progress(0, text="Procesando lote...")
                        total = len(df)
                        
                        for idx, row in df.iterrows():
                            c_name = str(row["name"]).strip()
                            f_group = str(row["functional_group"]).strip()
                            progress_bar.progress((idx + 1) / total, text=f"Procesando: {c_name}")
                            
                            smiles = name_to_smiles_simple(c_name)
                            if smiles:
                                img_b64 = smiles_to_png_base64(smiles)
                                st.session_state.dg_questions.append({
                                    "name": c_name,
                                    "functional_group": f_group,
                                    "smiles": smiles,
                                    "image_b64": img_b64
                                })
                                success_count += 1
                            else:
                                fail_count += 1
                                
                        progress_bar.empty()
                        st.success(t["bulk_success"].format(success_count, fail_count))
                except Exception as e:
                    st.error(f"Error procesando el archivo: {e}")

    # ================= COLUMNA DERECHA: LISTADO Y EXPORTACIÓN =================
    with col_list:
        st.subheader(t["section_export"])
        
        if st.session_state.dg_questions:
            st.write(f"Total de preguntas en la lista: **{len(st.session_state.dg_questions)}**")
            
            # --- BOTONES PRINCIPALES EN VERTICAL ---
            xml_data = generate_short_answer_xml(st.session_state.dg_questions, lang)
            
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
                file_name="moodle_functional_groups.xml",
                icon=":material/download:",
                mime="application/xml",
                use_container_width=True
            )

            if st.button(t["clear_all"], type="secondary", icon=":material/delete_sweep:", use_container_width=True):
                st.session_state.dg_questions = []
                st.session_state.last_fg_preview = None
                st.rerun()
            
            st.divider()
                
            # Contenedor con scroll para la lista de preguntas
            container_list = st.container(height=420)
            with container_list:
                for i, q in enumerate(st.session_state.dg_questions):
                    cols_item = st.columns([1, 2.5, 0.6])
                    with cols_item[0]:
                        if q['image_b64']:
                            st.image(f"data:image/png;base64,{q['image_b64']}", width=100)
                    with cols_item[1]:
                        st.markdown(f"**#{i+1}** (`{q['name']}`)<br>Grupo: **{q['functional_group']}**", unsafe_allow_html=True)
                    with cols_item[2]:
                        if st.button("🗑️", key=f"del_dg_{i}"):
                            st.session_state.dg_questions.pop(i)
                            st.rerun()
                    st.divider()
        else:
            st.info(t["empty_list"])