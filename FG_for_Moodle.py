# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:51:39 2026

@author: Carlos
"""

"""
Functional Group Determination Module for the Moodle Chemistry Suite.
Allows individual entry with automatic detection and manual editing (even if detection fails),
as well as bulk upload with detailed error reporting and cause tracking.
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

# Master dictionary with valid SMARTS strings and multiple accepted answers for ES and EN
FUNCTIONAL_GROUPS_MASTER = {
    'Ácido carboxílico': {
        'smarts': '[$([CX3](=O)[OX2H1]),$([CX3](=O)[O-])]', 
        'priority': 1,
        'es': ['Ácido carboxílico'],
        'en': ['Carboxylic acid']
    },
    'Peroxiácido (Ácido peroxicarboxílico)': {
        'smarts': '[CX3](=O)[OX2][OX2H]', 
        'priority': 2,
        'es': ['Peroxiácido (Ácido peroxicarboxílico)', 'Peroxiácido'],
        'en': ['Peroxyacid (Peroxycarboxylic acid)', 'Peroxyacid']
    },
    'Anhídrido de ácido': {
        'smarts': '[CX3](=O)[OX2][CX3](=O)', 
        'priority': 3,
        'es': ['Anhídrido de ácido', 'Anhídrido'],
        'en': ['Acid anhydride', 'Anhydride']
    },
    'Éster': {
        'smarts': '[#6][CX3](=O)[OX2][#6]', 
        'priority': 4,
        'es': ['Éster'],
        'en': ['Ester']
    },
    'Haluro de ácido (Acilo)': {
        'smarts': '[CX3](=O)[F,Cl,Br,I]', 
        'priority': 5,
        'es': ['Haluro de ácido (Acilo)', 'Haluro de acilo', 'Haluro de ácido'],
        'en': ['Acyl halide', 'Acid halide']
    },
    'Amida': {
        'smarts': '[CX3](=O)[NX3]', 
        'priority': 6,
        'es': ['Amida'],
        'en': ['Amide']
    },
    'Nitrilo': {
        'smarts': '[CX2,CX3]#[NX1]', 
        'priority': 7,
        'es': ['Nitrilo', 'Cianuro'],
        'en': ['Nitrile', 'Cyanide']
    },
    'Aldehído': {
        'smarts': '[CX3H1](=O)[#6]', 
        'priority': 8,
        'es': ['Aldehído'],
        'en': ['Aldehyde']
    },
    'Cetona': {
        'smarts': '[#6][CX3](=O)[#6]', 
        'priority': 9,
        'es': ['Cetona'],
        'en': ['Ketone']
    },
    'Imina / Base de Schiff': {
        'smarts': '[CX3](=[NX2][#6])[#6,H]', 
        'priority': 10,
        'es': ['Imina / Base de Schiff', 'Imina'],
        'en': ['Imine / Schiff base', 'Imine']
    },
    'Hemiacetal / Hemiketal': {
        'smarts': '[CX4]([OX2H])([OX2][#6])', 
        'priority': 11,
        'es': ['Hemiacetal / Hemiketal'],
        'en': ['Hemiacetal / Hemiketal']
    },
    'Alcohol': {
        'smarts': '[#6][OX2H]', 
        'priority': 12,
        'es': ['Alcohol', 'Grupo hidroxilo', 'Hidroxilo'],
        'en': ['Alcohol', 'Hydroxyl']
    },
    'Amina': {
        'smarts': '[NX3;H2,H1,H0;!$(NC=O)]', 
        'priority': 13,
        'es': ['Amina'],
        'en': ['Amine']
    },
    'Acetal / Ketal': {
        'smarts': '[CX4]([OX2][#6])([OX2][#6])', 
        'priority': 14,
        'es': ['Acetal / Ketal'],
        'en': ['Acetal / Ketal']
    },
    'Epóxido (Oxirano)': {
        'smarts': '[CX4]1[OX2][CX4]1', 
        'priority': 15,
        'es': ['Epóxido (Oxirano)', 'Epóxido', 'Oxirano'],
        'en': ['Epoxide (Oxirane)', 'Epoxide', 'Oxirane']
    },
    'Peróxido': {
        'smarts': '[OX2][OX2]', 
        'priority': 16,
        'es': ['Peróxido'],
        'en': ['Peroxide']
    },
    'Éter': {
        'smarts': '[#6][OX2][#6]', 
        'priority': 17,
        'es': ['Éter'],
        'en': ['Ether']
    },
    'Grupo Nitro': {
        'smarts': '[$([NX3](=O)=O),$([NX3+]([O-])=O)]', 
        'priority': 18,
        'es': ['Grupo Nitro', 'Nitro'],
        'en': ['Nitro group', 'Nitro']
    },
    'Haluro de alquilo': {
        'smarts': '[CX4][F,Cl,Br,I]', 
        'priority': 19,
        'es': ['Haluro de alquilo', 'Haloalcano'],
        'en': ['Alkyl halide', 'Haloalkane']
    },
    'Alquino (Triple enlace)': {
        'smarts': '[CX2]#[CX2]', 
        'priority': 20,
        'es': ['Alquino (Triple enlace)', 'Alquino', 'Triple enlace'],
        'en': ['Alkyne (Triple bond)', 'Alkyne', 'Triple bond']
    },
    'Alqueno (Doble enlace)': {
        'smarts': '[CX3]=[CX3]', 
        'priority': 21,
        'es': ['Alqueno (Doble enlace)', 'Alqueno', 'Doble enlace'],
        'en': ['Alkene (Double bond)', 'Alkene', 'Double bond']
    },
    'Anillo aromático (Bencénico)': {
        'smarts': 'c1ccccc1', 
        'priority': 22,
        'es': ['Anillo aromático (Bencénico)', 'Benceno', 'Anillo aromático', 'Aromático'],
        'en': ['Aromatic ring (Benzene)', 'Benzene', 'Aromatic ring']
    },
    'Alcano': {
        'smarts': '[CX4]', 
        'priority': 23,
        'es': ['Alcano'],
        'en': ['Alkane']
    }
}

# Bilingual interface texts
TEXTS_DG = {
    "es": {
        "title": "Determinación Automática de Grupo Funcional 🧪",
        "desc": "Introduce una molécula y RDKit detectará automáticamente su grupo funcional (incluyendo alcanos). Si lo deseas, puedes modificarlo antes de generar la pregunta.",
        "tab_individual": "Entrada Individual",
        "tab_bulk": "Carga Masiva",
        "section_export": "Listado y Exportación de Preguntas",
        "compound_name": "Nombre de la molécula (en inglés, ej: hexane, acetic acid)",
        "detect_btn": "Analizar y Detectar Grupo",
        "add_question": "Añadir a la Lista de Preguntas",
        "fg_label": "Grupo funcional (editable, separa sinónimos por comas si procede)",
        "bulk_uploader": "Sube un archivo Excel (.xlsx) o CSV",
        "bulk_note": "El archivo puede contener la columna **'name'** y opcionalmente **'functional_group'**. Si está vacía, se detectará automáticamente.",
        "start_bulk": "Procesar Carga Masiva",
        "question_title": "Identificación de Grupo Funcional: {}",
        "question_text": "¿Cuál es el grupo funcional principal de esta molécula?",
        "export_xml": "Descargar XML Moodle (Respuesta Corta)",
        "clear_all": "Borrar Toda la Lista",
        "error_rdkit": "RDKit no está disponible en este entorno para procesar las moléculas.",
        "error_smiles": "No se pudo encontrar un SMILES válido para el compuesto '{}' en NCI CIR.",
        "error_enter_name": "Por favor, introduce el nombre de la molécula.",
        "error_fg_empty": "El grupo funcional no puede estar vacío.",
        "error_col_name": "El archivo debe contener al menos una columna llamada 'name' o 'nombre'.",
        "success_analyzed": "Molécula analizada: `{}`",
        "spinner_analyzing": "Buscando estructura y analizando grupo funcional...",
        "progress_batch": "Procesando lote y detectando grupos...",
        "success_added": "Pregunta añadida correctamente con el grupo: **{}**.",
        "bulk_success": "Carga masiva completada: **{}** éxitos, **{}** fallos.",
        "bulk_errors_title": "Ver detalle de errores en la carga masiva",
        "empty_list": "Aún no hay preguntas añadidas. Completa el formulario de la izquierda.",
        "last_added": "Última molécula añadida",
        "group_label": "Grupo",
        "group_detected": "Grupo asignado"
    },
    "en": {
        "title": "Automatic Functional Group Determination 🧪",
        "desc": "Enter a molecule and RDKit will automatically detect its functional group (including alkanes). You can edit it before generating the question.",
        "tab_individual": "Individual Entry",
        "tab_bulk": "Bulk Upload",
        "section_export": "Questions List & Export",
        "compound_name": "Molecule Name (in English, e.g., hexane, acetic acid)",
        "detect_btn": "Analyze and Detect Group",
        "add_question": "Add to Questions List",
        "fg_label": "Functional group (editable)",
        "bulk_uploader": "Upload Excel (.xlsx) or CSV file",
        "bulk_note": "The file can contain a **'name'** column and optionally **'functional_group'**. If empty, automatic detection will be applied.",
        "start_bulk": "Process Bulk Upload",
        "question_title": "Functional Group Identification: {}",
        "question_text": "What is the main functional group of this molecule?",
        "export_xml": "Download Moodle XML (Short Answer)",
        "clear_all": "Clear All List",
        "error_rdkit": "RDKit is not available in this environment to process molecules.",
        "error_smiles": "Could not find a valid SMILES for compound '{}' via NCI CIR.",
        "error_enter_name": "Please enter a molecule name.",
        "error_fg_empty": "The functional group cannot be empty.",
        "error_col_name": "The file must contain at least a 'name' column.",
        "success_analyzed": "Analyzed molecule: `{}`",
        "spinner_analyzing": "Searching structure and analyzing functional group...",
        "progress_batch": "Processing batch and detecting groups...",
        "success_added": "Question successfully added with group: **{}**.",
        "bulk_success": "Bulk processing complete: **{}** successes, **{}** failures.",
        "bulk_errors_title": "View bulk upload error details",
        "empty_list": "No questions added yet. Fill out the form on the left.",
        "last_added": "Last added molecule",
        "group_label": "Group",
        "group_detected": "Assigned group"
    }
}

def detectar_grupo_funcional(mol, lang="es"):
    """
    Detects the main functional group using the master dictionary sorted by IUPAC priority.
    Returns the full list of accepted answers/synonyms for the active language.
    If no specific group matches, defaults to Alkane.
    """
    if not mol:
        return ['Alcano'] if lang == 'es' else ['Alkane']

    grupos_ordenados = sorted(FUNCTIONAL_GROUPS_MASTER.items(), key=lambda x: x[1]['priority'])

    for nombre_es, data in grupos_ordenados:
        patron = Chem.MolFromSmarts(data['smarts'])
        if patron and mol.HasSubstructMatch(patron):
            return data['es'] if lang == 'es' else data['en']

    # Fallback to Alkane
    return ['Alcano'] if lang == 'es' else ['Alkane']

def name_to_smiles_simple(compound_name):
    """Searches for a molecule's SMILES using the NCI CIR API."""
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
                        return Chem.MolToSmiles(mol, canonical=True), mol
                return smiles, None
    except Exception:
        pass
    return None, None

def smiles_to_png_base64(mol):
    """Generates a Base64 PNG image of the molecule using clean RDKit rendering."""
    if not RDKIT_AVAILABLE or not mol:
        return None
    try:
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
    """Generates Moodle XML for short-answer questions including base64-embedded images."""
    quiz = ET.Element("quiz")
    
    for item in questions_list:
        name = item["name"]
        fg = item["functional_group"]
        b64_img = item["image_b64"]
        
        question = ET.Element("question", type="shortanswer")
        
        q_name = ET.SubElement(question, "name")
        q_text_name = ET.SubElement(q_name, "text")
        q_text_name.text = TEXTS_DG[lang]["question_title"].format(name)
        
        q_text = ET.SubElement(question, "questiontext", format="html")
        q_text_content = ET.SubElement(q_text, "text")
        
        base_prompt = TEXTS_DG[lang]["question_text"]
        img_tag = f'<p><img src="data:image/png;base64,{b64_img}" alt="{name}" /></p>' if b64_img else ""
        q_text_content.text = f"<p>{base_prompt}</p>{img_tag}"
        
        # Handle functional groups (supports both single string or list of valid answers)
        answers_list = fg if isinstance(fg, list) else [fg]
        
        for ans_val in answers_list:
            answer = ET.SubElement(question, "answer", fraction="100", format="moodle_auto_format")
            ans_text = ET.SubElement(answer, "text")
            ans_text.text = ans_val.strip()
            
            feedback = ET.SubElement(answer, "feedback", format="html")
            ET.SubElement(feedback, "text").text = "¡Correcto!" if lang == "es" else "Correct!"
        
        quiz.append(question)
        
    tree = ET.ElementTree(quiz)
    if hasattr(ET, "indent"):
        ET.indent(tree, space="    ")
    return ET.tostring(tree.getroot(), encoding="utf-8", xml_declaration=True)

def render_nueva_herramienta(texts_main):
    """Main function that renders the module in Streamlit."""
    lang = st.session_state.get("lang", "en")
    t = TEXTS_DG[lang]
    
    st.title(t["title"])
    st.write(t["desc"])
    
    if not RDKIT_AVAILABLE:
        st.warning(t["error_rdkit"])
        return

    if "dg_questions" not in st.session_state:
        st.session_state.dg_questions = []
        
    if "last_fg_preview" not in st.session_state:
        st.session_state.last_fg_preview = None

    if "temp_fg_data" not in st.session_state:
        st.session_state.temp_fg_data = None

    # --- TWO-COLUMN LAYOUT ---
    col_input, col_list = st.columns([1.1, 0.9], gap="large")

    # ================= LEFT COLUMN: INPUT (INDIVIDUAL & BULK) =================
    with col_input:
        tab_ind, tab_bulk = st.tabs([t["tab_individual"], t["tab_bulk"]])

        with tab_ind:
            # Step 1: Input molecule and analyze
            with st.form("form_individual_dg_detect"):
                c_name = st.text_input(t["compound_name"])
                submitted_detect = st.form_submit_button(t["detect_btn"], type="primary", icon=":material/auto_fix_high:")
                
                if submitted_detect:
                    if not c_name.strip():
                        st.error(t["error_enter_name"])
                        st.session_state.temp_fg_data = None
                    else:
                        with st.spinner(t["spinner_analyzing"]):
                            smiles, mol = name_to_smiles_simple(c_name)
                            if not smiles or not mol:
                                st.error(t["error_smiles"].format(c_name))
                                st.session_state.temp_fg_data = None
                            else:
                                detected_fg_list = detectar_grupo_funcional(mol, lang)
                                
                                img_b64 = smiles_to_png_base64(mol)
                                st.session_state.temp_fg_data = {
                                    "name": c_name.strip(),
                                    "functional_group": detected_fg_list, # Es una lista de sinónimos
                                    "smiles": smiles,
                                    "image_b64": img_b64
                                }

            # Step 2: If analyzed, show preview, allow manual editing and confirm
            if st.session_state.temp_fg_data:
                data = st.session_state.temp_fg_data
                st.success(t["success_analyzed"].format(data['name']))
                
                if data["image_b64"]:
                    st.image(f"data:image/png;base64,{data['image_b64']}", width=250, caption=data["name"])

                with st.form("form_individual_dg_confirm"):
                    # Mostramos los sinónimos unidos por coma en el input editable por comodidad
                    current_fg_str = ", ".join(data["functional_group"]) if isinstance(data["functional_group"], list) else data["functional_group"]
                    edited_fg = st.text_input(t["fg_label"], value=current_fg_str)
                    submitted_add = st.form_submit_button(t["add_question"], type="primary", icon=":material/add_task:")
                    
                    if submitted_add:
                        if not edited_fg.strip():
                            st.error(t["error_fg_empty"])
                        else:
                            # Convertimos lo introducido de nuevo en lista separando por comas
                            final_fg_list = [x.strip() for x in edited_fg.split(",") if x.strip()]
                            
                            st.session_state.dg_questions.append({
                                "name": data["name"],
                                "functional_group": final_fg_list,
                                "smiles": data["smiles"],
                                "image_b64": data["image_b64"]
                            })
                            st.session_state.last_fg_preview = {
                                "name": data["name"],
                                "functional_group": " / ".join(final_fg_list),
                                "image_b64": data["image_b64"]
                            }
                            st.session_state.temp_fg_data = None
                            st.success(t["success_added"].format(", ".join(final_fg_list)))
                            st.rerun()

            if st.session_state.last_fg_preview:
                prev = st.session_state.last_fg_preview
                st.markdown(f"**{t['last_added']}:** `{prev['name']}`<br>{t['group_detected']}: **{prev['functional_group']}**", unsafe_allow_html=True)

        with tab_bulk:
            st.info(t["bulk_note"])
            uploaded_file = st.file_uploader(t["bulk_uploader"], type=["xlsx", "csv"], key="bulk_dg_file")
            
            if uploaded_file and st.button(t["start_bulk"]):
                try:
                    if uploaded_file.name.endswith(".csv"):
                        df = pd.read_csv(uploaded_file)
                    else:
                        df = pd.read_excel(uploaded_file)
                        
                    col_name = 'name' if 'name' in df.columns else ('nombre' if 'nombre' in df.columns else None)
                    if not col_name:
                        st.error(t["error_col_name"])
                    else:
                        col_fg = 'functional_group' if 'functional_group' in df.columns else ('grupo_funcional' if 'grupo_funcional' in df.columns else None)
                        
                        success_count = 0
                        fail_count = 0
                        errores_detalle = []
                        
                        progress_bar = st.progress(0, text=t["progress_batch"])
                        total = len(df)
                        
                        for idx, row in df.iterrows():
                            c_name = str(row[col_name]).strip()
                            progress_bar.progress((idx + 1) / total, text=f"Processing: {c_name}")
                            
                            custom_fg_list = None
                            if col_fg and pd.notna(row[col_fg]) and str(row[col_fg]).strip():
                                raw_custom = str(row[col_fg]).strip()
                                custom_fg_list = [x.strip() for x in raw_custom.split(",") if x.strip()]
                            
                            smiles, mol = name_to_smiles_simple(c_name)
                            if not smiles or not mol:
                                fail_count += 1
                                errores_detalle.append({"row": idx + 1, "name": c_name, "cause": "SMILES not found in NCI CIR"})
                                continue
                                
                            if custom_fg_list:
                                f_group = custom_fg_list
                            else:
                                f_group = detectar_grupo_funcional(mol, lang)
                                
                            img_b64 = smiles_to_png_base64(mol)
                            st.session_state.dg_questions.append({
                                "name": c_name,
                                "functional_group": f_group,
                                "smiles": smiles,
                                "image_b64": img_b64
                            })
                            success_count += 1
                                
                        progress_bar.empty()
                        st.success(t["bulk_success"].format(success_count, fail_count))
                        
                        if errores_detalle:
                            with st.expander(t["bulk_errors_title"]):
                                df_err = pd.DataFrame(errores_detalle)
                                st.dataframe(df_err, use_container_width=True)
                                
                except Exception as e:
                    st.error(f"Error processing file: {e}")

    # ================= COLUMNA DERECHA: LISTADO Y EXPORTACIÓN =================
    with col_list:
        st.subheader(t["section_export"])
        
        if st.session_state.dg_questions:
            st.write(f"Total questions in list: **{len(st.session_state.dg_questions)}**")
            
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
                st.session_state.temp_fg_data = None
                st.rerun()
            
            st.divider()
                
            container_list = st.container(height=420)
            with container_list:
                for i, q in enumerate(st.session_state.dg_questions):
                    cols_item = st.columns([1, 2.5, 0.6])
                    with cols_item[0]:
                        if q['image_b64']:
                            st.image(f"data:image/png;base64,{q['image_b64']}", width=100)
                    with cols_item[1]:
                        fg_display = " / ".join(q['functional_group']) if isinstance(q['functional_group'], list) else q['functional_group']
                        st.markdown(f"**#{i+1}** (`{q['name']}`)<br>{t['group_label']}: **{fg_display}**", unsafe_allow_html=True)
                    with cols_item[2]:
                        if st.button("🗑️", key=f"del_dg_{i}"):
                            st.session_state.dg_questions.pop(i)
                            st.rerun()
                    st.divider()
        else:
            st.info(t["empty_list"])
