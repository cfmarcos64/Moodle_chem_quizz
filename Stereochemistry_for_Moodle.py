# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 10:15:20 2026

@author: Carlos Fernandez
"""

"""
Stereochemistry and Chiral/Geometric Centers Module for the Moodle Chemistry Suite.
Allows automatic detection of chiral centers (R/S) and double bonds (Z/E) with RDKit,
adapted prompt texts, two-column layout, live preview, XML export, and detailed bulk upload error reporting.
"""

import streamlit as str_lit
import xml.etree.ElementTree as ET
import pandas as pd
import requests
import io
import base64

# --- Dependency Verification (RDKit) ---
try:
    from rdkit import Chem
    from rdkit.Chem import FindMolChiralCenters, Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Bilingual interface texts and adapted prompts
TEXTS_EQ = {
    "es": {
        "title": "Generador de Preguntas de Estereoquímica 🧬",
        "desc": "Detecta automáticamente centros quirales (R/S) y dobles enlaces (Z/E) con RDKit, indica sus índices en la molécula y genera preguntas de respuesta corta para Moodle.",
        "tab_individual": "Entrada Individual",
        "tab_bulk": "Carga Masiva",
        "section_export": "Listado y Exportación de Preguntas",
        "compound_name": "Nombre de la molécula en inglés (ej: 2-butene, L-lactic acid)",
        "add_question": "Analizar y Añadir a la Lista",
        "bulk_uploader": "Sube un archivo Excel (.xlsx) o CSV",
        "bulk_note": "El archivo debe contener al menos una columna llamada **'name'** o **'nombre'** con los nombres en inglés.",
        "start_bulk": "Procesar Carga Masiva",
        "question_title": "Estereoquímica de {}",
        "prompt_chiral": "Indica la configuración absoluta de los centros quirales señalados (formato: índiceR/S, ej: 1S o 2R, 4S):",
        "prompt_geometric": "Indica la configuración geométrica de los dobles enlaces señalados (formato: índiceZ/E, ej: 2Z o 3E):",
        "prompt_both": "Indica la configuración estereoquímica absoluta y geométrica de los elementos señalados (ej: 2R, 4E):",
        "export_xml": "Descargar XML Moodle (Respuesta Corta)",
        "clear_all": "Borrar Toda la Lista",
        "error_rdkit": "RDKit no está disponible en este entorno para realizar los cálculos estereoquímicos.",
        "error_smiles": "No se pudo encontrar un SMILES válido para el compuesto '{}' en NCI CIR.",
        "error_chiral": "La molécula analizada no presenta elementos estereogénicos detectables (centros quirales ni dobles enlaces Z/E).",
        "error_enter_name": "Por favor, introduce el nombre de la molécula.",
        "error_col_name": "El archivo debe contener al menos una columna llamada 'name' o 'nombre'.",
        "success_added": "Pregunta de estereoquímica añadida correctamente.",
        "bulk_success": "Carga masiva completada: **{}** éxitos, **{}** fallos.",
        "bulk_errors_title": "Ver detalle de errores en la carga masiva",
        "empty_list": "Aún no hay preguntas añadidas. Completa el formulario de la izquierda.",
        "last_added": "Última molécula añadida",
        "config_label": "Config"
    },
    "en": {
        "title": "Stereochemistry Question Generator 🧬",
        "desc": "Automatically detects chiral centers (R/S) and double bonds (Z/E) with RDKit, indicates their indices on the molecule, and generates short-answer Moodle questions.",
        "tab_individual": "Individual Entry",
        "tab_bulk": "Bulk Upload",
        "section_export": "Questions List & Export",
        "compound_name": "Molecule name in English (e.g., 2-butene, L-lactic acid)",
        "add_question": "Analyze and Add to List",
        "bulk_uploader": "Upload Excel (.xlsx) or CSV file",
        "bulk_note": "The file must contain a column named **'name'** or **'nombre'** with molecule names in English.",
        "start_bulk": "Process Bulk Upload",
        "question_title": "Stereochemistry of {}",
        "prompt_chiral": "Indicate the absolute configuration of the specified chiral centers (format: indexR/S, e.g., 1S or 2R, 4S):",
        "prompt_geometric": "Indicate the geometric configuration of the specified double bonds (format: indexZ/E, e.g., 2Z or 3E):",
        "prompt_both": "Indicate the absolute and geometric stereochemical configuration of the specified elements (e.g., 2R, 4E):",
        "export_xml": "Download Moodle XML (Short Answer)",
        "clear_all": "Clear All Questions",
        "error_rdkit": "RDKit is not available in this environment to perform stereochemical calculations.",
        "error_smiles": "Could not find a valid SMILES for compound '{}' via NCI CIR.",
        "error_chiral": "The analyzed molecule has no detectable stereogenic elements (chiral centers or Z/E double bonds).",
        "error_enter_name": "Please enter a molecule name.",
        "error_col_name": "The file must contain at least a 'name' or 'nombre' column.",
        "success_added": "Stereochemistry question successfully added.",
        "bulk_success": "Bulk processing complete: **{}** successes, **{}** failures.",
        "bulk_errors_title": "View bulk upload error details",
        "empty_list": "No questions added yet. Fill out the form on the left.",
        "last_added": "Last added molecule",
        "config_label": "Config"
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

def analizar_estereoquimica_mol(nombre_molecula):
    """
    Loads the molecule, calculates CIP stereochemistry (chiral and Z/E double bonds),
    assigns labels/notes to corresponding atoms and bonds, and generates the image.
    """
    if not RDKIT_AVAILABLE:
        return None, None, [], "RDKit not available"

    smiles = buscar_smiles_nci(nombre_molecula)
    if not smiles:
        return None, None, [], f"SMILES not found for {nombre_molecula}"

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None, None, [], f"Error parsing SMILES: {smiles}"

    try:
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    except Exception as e:
        return None, None, [], f"Error assigning stereochemistry: {e}"

    elementos_estereo = []
    indices_resaltar_atomos = []
    indices_resaltar_enlaces = []

    # 1. Search for chiral centers (R/S)
    centros_quirales = FindMolChiralCenters(mol, includeUnassigned=True)
    for idx, config in centros_quirales:
        if config in ["R", "S"]:
            elementos_estereo.append((idx, f"{idx}{config}", "atom"))
            indices_resaltar_atomos.append(idx)
            atom = mol.GetAtomWithIdx(idx)
            atom.SetProp('atomNote', f"{idx}")

    # 2. Search for double bonds with geometric isomerism (Z/E)
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            stereo = bond.GetStereo()
            if stereo in [Chem.rdchem.BondStereo.STEREOZ, Chem.rdchem.BondStereo.STEREOE]:
                label_ze = "Z" if stereo == Chem.rdchem.BondStereo.STEREOZ else "E"
                begin_idx = bond.GetBeginAtomIdx()
                b_idx = bond.GetIdx()
                elementos_estereo.append((begin_idx, f"{begin_idx}{label_ze}", "bond"))
                indices_resaltar_enlaces.append(b_idx)
                begin_atom = mol.GetAtomWithIdx(begin_idx)
                begin_atom.SetProp('atomNote', f"{begin_idx}")

    if not elementos_estereo:
        return None, None, [], "No detectable stereogenic elements"

    elementos_estereo.sort(key=lambda x: x[0])
    respuestas_formateadas = [item[1] for item in elementos_estereo]
    respuesta_correcta = ", ".join(respuestas_formateadas)

    has_chiral = any(t == "atom" for _, _, t in elementos_estereo)
    has_geometric = any(t == "bond" for _, _, t in elementos_estereo)
    
    if has_chiral and has_geometric:
        tipo_estereo = "both"
    elif has_geometric:
        tipo_estereo = "geometric"
    else:
        tipo_estereo = "chiral"

    opts = Draw.MolDrawOptions()
    opts.addAtomIndices = False  
    opts.fixedFontSize = 14

    try:
        img = Draw.MolToImage(
            mol,
            size=(500, 400),
            highlightAtoms=indices_resaltar_atomos,
            highlightBonds=indices_resaltar_enlaces,
            drawOptions=opts
        )
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
        img_b64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
    except Exception:
        img_b64 = None

    return respuesta_correcta, img_b64, tipo_estereo, None

def generate_chiral_short_answer_xml(questions_list, lang):
    """Generates Moodle XML for short-answer questions adapting the prompt based on stereo type."""
    quiz = ET.Element("quiz")
    
    for item in questions_list:
        name = item["name"]
        correct_ans = item["correct_answer"]
        b64_img = item["image_b64"]
        tipo_estereo = item.get("tipo_estereo", "chiral")
        
        question = ET.Element("question", type="shortanswer")
        
        q_name = ET.SubElement(question, "name")
        q_text_name = ET.SubElement(q_name, "text")
        q_text_name.text = TEXTS_EQ[lang]["question_title"].format(name)
        
        if tipo_estereo == "geometric":
            base_prompt = TEXTS_EQ[lang]["prompt_geometric"]
        elif tipo_estereo == "both":
            base_prompt = TEXTS_EQ[lang]["prompt_both"]
        else:
            base_prompt = TEXTS_EQ[lang]["prompt_chiral"]
        
        q_text = ET.SubElement(question, "questiontext", format="html")
        q_text_content = ET.SubElement(q_text, "text")
        
        img_tag = f'<p><img src="data:image/png;base64,{b64_img}" alt="{name}" /></p>' if b64_img else ""
        q_text_content.text = f"<![CDATA[<p>{base_prompt}</p> {img_tag}]]>"
        
        answer = ET.SubElement(question, "answer", fraction="100", format="moodle_auto_format")
        ans_text = ET.SubElement(answer, "text")
        ans_text.text = correct_ans
        
        feedback = ET.SubElement(answer, "feedback", format="html")
        ET.SubElement(feedback, "text").text = "¡Correcto!" if lang == "es" else "Correct!"
        
        quiz.append(question)
        
    tree = ET.ElementTree(quiz)
    ET.indent(tree, space="  ")
    return ET.tostring(tree.getroot(), encoding="utf-8", xml_declaration=True)

def render_nueva_herramienta(texts_main):
    """Main function that renders the stereochemistry module in Streamlit."""
    lang = str_lit.session_state.get("lang", "en")
    t = TEXTS_EQ[lang]
    
    str_lit.title(t["title"])
    str_lit.write(t["desc"])
    
    if not RDKIT_AVAILABLE:
        str_lit.warning(t["error_rdkit"])
        return

    if "eq_questions" not in str_lit.session_state:
        str_lit.session_state.eq_questions = []
        
    if "last_eq_preview" not in str_lit.session_state:
        str_lit.session_state.last_eq_preview = None

    # --- TWO-COLUMN LAYOUT ---
    col_input, col_list = str_lit.columns([1.1, 0.9], gap="large")

    # ================= LEFT COLUMN: INPUT (INDIVIDUAL & BULK) =================
    with col_input:
        tab_ind, tab_bulk = str_lit.tabs([t["tab_individual"], t["tab_bulk"]])

        with tab_ind:
            with str_lit.form("form_individual_eq", clear_on_submit=True):
                c_name = str_lit.text_input(t["compound_name"])
                submitted = str_lit.form_submit_button(t["add_question"], type="primary", icon=":material/add_task:")
                
                if submitted:
                    if not c_name.strip():
                        str_lit.error(t["error_enter_name"])
                    else:
                        with str_lit.spinner("Searching structure, calculating stereochemistry (R/S & Z/E), and rendering image..." if lang=="en" else "Buscando estructura, calculando estereoquímica (R/S y Z/E) y renderizando imagen..."):
                            correct_ans, img_b64, tipo_estereo, err = analizar_estereoquimica_mol(c_name)
                            if err:
                                str_lit.error(f"Error: {err}")
                            elif not correct_ans:
                                str_lit.warning(t["error_chiral"])
                            else:
                                str_lit.session_state.eq_questions.append({
                                    "name": c_name.strip(),
                                    "correct_answer": correct_ans,
                                    "image_b64": img_b64,
                                    "tipo_estereo": tipo_estereo
                                })
                                str_lit.session_state.last_eq_preview = {
                                    "name": c_name.strip(),
                                    "image_b64": img_b64,
                                    "correct_answer": correct_ans
                                }
                                str_lit.success(f"{t['success_added']} Answer: **{correct_ans}**" if lang=="en" else f"{t['success_added']} Respuesta: **{correct_ans}**")

            if str_lit.session_state.last_eq_preview:
                prev = str_lit.session_state.last_eq_preview
                str_lit.markdown(f"**{t['last_added']}:** `{prev['name']}` (Config: `{prev['correct_answer']}`)")
                if prev["image_b64"]:
                    str_lit.image(f"data:image/png;base64,{prev['image_b64']}", width=250, caption=prev["name"])

        with tab_bulk:
            str_lit.info(t["bulk_note"])
            uploaded_file = str_lit.file_uploader(t["bulk_uploader"], type=["xlsx", "csv"], key="uploader_eq")
            
            if uploaded_file and str_lit.button(t["start_bulk"]):
                try:
                    if uploaded_file.name.endswith(".csv"):
                        df = pd.read_csv(uploaded_file)
                    else:
                        df = pd.read_excel(uploaded_file)
                        
                    col_name = 'name' if 'name' in df.columns else ('nombre' if 'nombre' in df.columns else None)
                    
                    if not col_name:
                        str_lit.error(t["error_col_name"])
                    else:
                        success_count = 0
                        fail_count = 0
                        errores_detalle = []
                        
                        progress_bar = str_lit.progress(0, text="Processing batch..." if lang=="en" else "Procesando lote estereoquímico...")
                        total = len(df)
                        
                        for idx, row in df.iterrows():
                            c_name = str(row[col_name]).strip() if pd.notna(row[col_name]) else ""
                            progress_bar.progress((idx + 1) / total, text=f"Analyzing: {c_name}")
                            
                            if not c_name or c_name == "nan":
                                fail_count += 1
                                errores_detalle.append({"row": idx + 1, "name": "N/A", "cause": "Missing molecule name"})
                                continue
                                
                            correct_ans, img_b64, tipo_estereo, err = analizar_estereoquimica_mol(c_name)
                            if err:
                                fail_count += 1
                                errores_detalle.append({"row": idx + 1, "name": c_name, "cause": err})
                            elif not correct_ans:
                                fail_count += 1
                                errores_detalle.append({"row": idx + 1, "name": c_name, "cause": "No detectable stereogenic elements"})
                            else:
                                str_lit.session_state.eq_questions.append({
                                    "name": c_name,
                                    "correct_answer": correct_ans,
                                    "image_b64": img_b64,
                                    "tipo_estereo": tipo_estereo
                                })
                                success_count += 1
                                
                        progress_bar.empty()
                        str_lit.success(t["bulk_success"].format(success_count, fail_count))
                        
                        if errores_detalle:
                            with str_lit.expander(t["bulk_errors_title"]):
                                df_err = pd.DataFrame(errores_detalle)
                                str_lit.dataframe(df_err, use_container_width=True)
                                
                except Exception as e:
                    str_lit.error(f"Error processing file: {e}")

    # ================= RIGHT COLUMN: LIST & EXPORT =================
    with col_list:
        str_lit.subheader(t["section_export"])
        
        if str_lit.session_state.eq_questions:
            str_lit.write(f"Total questions in list: **{len(str_lit.session_state.eq_questions)}**" if lang=="en" else f"Total de preguntas en la lista: **{len(str_lit.session_state.eq_questions)}**")
            
            # --- VERTICALLY ALIGNED MAIN BUTTONS ---
            xml_data = generate_chiral_short_answer_xml(str_lit.session_state.eq_questions, lang)
            
            str_lit.markdown("""
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
            
            str_lit.download_button(
                label=t["export_xml"],
                data=xml_data,
                file_name="moodle_stereochemistry_questions.xml",
                icon=":material/download:",
                mime="application/xml",
                use_container_width=True
            )

            if str_lit.button(t["clear_all"], type="secondary", icon=":material/delete_sweep:", use_container_width=True):
                str_lit.session_state.eq_questions = []
                str_lit.session_state.last_eq_preview = None
                str_lit.rerun()
            
            str_lit.divider()
                
            # Scrollable container for questions list
            container_list = str_lit.container(height=420)
            with container_list:
                for i, q in enumerate(str_lit.session_state.eq_questions):
                    cols_item = str_lit.columns([1, 2.5, 0.6])
                    with cols_item[0]:
                        if q['image_b64']:
                            str_lit.image(f"data:image/png;base64,{q['image_b64']}", width=100)
                    with cols_item[1]:
                        str_lit.markdown(f"**#{i+1}** (`{q['name']}`)<br>{t['config_label']}: **`{q['correct_answer']}`**", unsafe_allow_html=True)
                    with cols_item[2]:
                        if str_lit.button("🗑️", key=f"del_eq_{i}"):
                            str_lit.session_state.eq_questions.pop(i)
                            str_lit.rerun()
                    str_lit.divider()
        else:
            str_lit.info(t["empty_list"])