# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 09:05:18 2026

@author: Carlos Fernandez
"""

"""
Módulo de Análisis de Estereoisómeros para la Suite de Moodle.
Permite comparar la molécula base (con centros quirales y/o alquenos Z/E) con estereoisómeros
automáticos o con otro compuesto introducido manualmente, clasificándolo en las 4 categorías requeridas.
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
    from rdkit.Chem import FindMolChiralCenters, Draw
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Textos bilingües para la interfaz
TEXTS_STEREO = {
    "es": {
        "title": "Análisis y Relación de Estereoisómeros 🔄",
        "desc": "Genera preguntas Moodle de opción múltiple comparando una molécula base (quiral o con alquenos Z/E) con sus estereoisómeros o con otra molécula manual.",
        "tab_individual": "Entrada Individual",
        "tab_bulk": "Carga Masiva",
        "section_export": "Listado y Exportación de Preguntas",
        "compound_name": "Nombre de la molécula base (en inglés, ej: 2-butene, D-glucose)",
        "btn_analyze": "Analizar Molécula Base",
        "mode_selection": "Selecciona el origen de la segunda molécula a comparar:",
        "option_auto": "Elegir estereoisómero generado automáticamente (quiral / Z,E)",
        "option_manual": "Introducir otro compuesto manualmente (nombre o SMILES)",
        "select_isomer": "Selecciona el estereoisómero:",
        "manual_comparison_name": "Nombre o SMILES del segundo compuesto:",
        "btn_evaluate_manual": "Evaluar Comparación Manual",
        "relationship_label": "Relación estereoquímica / tipo detectado:",
        "add_question": "Añadir Pregunta a la Lista",
        "bulk_uploader": "Sube un archivo Excel (.xlsx) o CSV",
        "bulk_note": "El archivo debe contener las columnas: **'molecule1'** y **'molecule2'**.",
        "start_bulk": "Procesar Carga Masiva",
        "question_title": "Relación Estereoquímica: {} vs {}",
        "question_text": "¿Cuál es la relación entre la molécula 1 (izquierda) y la molécula 2 (derecha)?",
        "export_xml": "Descargar XML Moodle (Opción Múltiple)",
        "clear_all": "Borrar Toda la Lista",
        "error_rdkit": "RDKit no está disponible en este entorno.",
        "error_smiles": "No se pudo encontrar un SMILES válido para '{}' en NCI CIR.",
        "success_added": "Pregunta añadida correctamente a la lista.",
        "bulk_success": "Carga masiva completada: **{}** éxitos, **{}** fallos.",
        "empty_list": "Aún no hay preguntas añadidas. Completa el formulario de la izquierda.",
        "rel_enantiomer": "Son enantiómeros",
        "rel_diastereomer": "Son diastereoisómeros",
        "rel_identical": "Son compuestos idénticos",
        "rel_different": "Son compuestos diferentes no estereoisoméricos"
    },
    "en": {
        "title": "Stereoisomer Analysis & Relationship 🔄",
        "desc": "Generate Moodle multiple-choice questions comparing a base molecule (chiral or Z/E alkenes) with its stereoisomers or a manual compound.",
        "tab_individual": "Individual Entry",
        "tab_bulk": "Bulk Upload",
        "section_export": "Questions List & Export",
        "compound_name": "Base molecule name (in English, e.g., 2-butene, D-glucose)",
        "btn_analyze": "Analyze Base Molecule",
        "mode_selection": "Select the origin of the second molecule to compare:",
        "option_auto": "Choose automatically generated stereoisomer (chiral / Z,E)",
        "option_manual": "Enter another compound manually (name or SMILES)",
        "select_isomer": "Select stereoisomer:",
        "manual_comparison_name": "Name or SMILES of the second compound:",
        "btn_evaluate_manual": "Evaluate Manual Comparison",
        "relationship_label": "Detected relationship / type:",
        "add_question": "Add Question to List",
        "bulk_uploader": "Upload Excel (.xlsx) or CSV file",
        "bulk_note": "The file must contain columns: **'molecule1'** and **'molecule2'**.",
        "start_bulk": "Process Bulk Upload",
        "question_title": "Stereochemical Relationship: {} vs {}",
        "question_text": "What is the relationship between molecule 1 (left) and molecule 2 (right)?",
        "export_xml": "Download Moodle XML (Multiple Choice)",
        "clear_all": "Clear All List",
        "error_rdkit": "RDKit is not available in this environment.",
        "error_smiles": "Could not find a valid SMILES for '{}' via NCI CIR.",
        "success_added": "Question successfully added to the list.",
        "bulk_success": "Bulk processing complete: **{}** successes, **{}** failures.",
        "empty_list": "No questions added yet. Fill out the form on the left.",
        "rel_enantiomer": "They are enantiomers",
        "rel_diastereomer": "They are diastereomers",
        "rel_identical": "They are identical compounds",
        "rel_different": "They are different non-stereoisomeric compounds"
    }
}

def resolver_smiles(input_str):
    """Obtiene el SMILES canónico a partir del nombre (vía NCI CIR) o valida un SMILES directo."""
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

def generar_imagen_comparativa_b64(mol1, mol2):
    """Genera una imagen PNG en Base64 con las dos moléculas lado a lado de forma limpia."""
    if not RDKIT_AVAILABLE:
        return None
    try:
        opts = Draw.MolDrawOptions()
        opts.addAtomIndices = False
        
        img = Draw.MolsToGridImage(
            [mol1, mol2],
            molsPerRow=2,
            subImgSize=(350, 350),
            legends=None,
            drawOptions=opts
        )
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
        return base64.b64encode(buffered.getvalue()).decode("utf-8")
    except Exception:
        pass
    return None

def obtener_descripcion_estereo(mol):
    """Obtiene una cadena descriptiva de los centros quirales y dobles enlaces Z/E de la molécula."""
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    partes = []
    
    chiral_centers = FindMolChiralCenters(mol, includeUnassigned=True)
    if chiral_centers:
        partes.append(", ".join([f"C{idx}:{conf}" for idx, conf in chiral_centers]))
        
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            stereo = bond.GetStereo()
            if stereo in [Chem.rdchem.BondStereo.STEREOZ, Chem.rdchem.BondStereo.STEREOE]:
                label = "Z" if stereo == Chem.rdchem.BondStereo.STEREOZ else "E"
                partes.append(f"Enlace {bond.GetBeginAtomIdx()}={bond.GetEndAtomIdx()}:{label}")
                
    return " | ".join(partes) if partes else "Sin estereocentros explícitos"

def calcular_relacion_estereoquimica(mol_orig, mol_sel):
    """
    Evalúa la relación exacta entre dos moléculas:
    - 'different': Si no tienen la misma fórmula o conectividad.
    - 'identical': Si son idénticas o forman una estructura meso.
    - 'enantiomer': Si son enantiómeros.
    - 'diastereomer': Si son diastereoisómeros (incluye isómeros Z/E).
    """
    try:
        smi_base_no_stereo = Chem.MolToSmiles(mol_orig, isomericSmiles=False)
        smi_comp_no_stereo = Chem.MolToSmiles(mol_sel, isomericSmiles=False)
        
        if smi_base_no_stereo != smi_comp_no_stereo:
            return "different"

        Chem.AssignStereochemistry(mol_orig, force=True, cleanIt=True)
        Chem.AssignStereochemistry(mol_sel, force=True, cleanIt=True)
        
        smi_orig_iso = Chem.MolToSmiles(mol_orig, isomericSmiles=True)
        smi_sel_iso = Chem.MolToSmiles(mol_sel, isomericSmiles=True)
        
        if smi_orig_iso == smi_sel_iso:
            return "identical"
            
        orig_dict = dict(FindMolChiralCenters(mol_orig, includeUnassigned=True))
        sel_dict = dict(FindMolChiralCenters(mol_sel, includeUnassigned=True))
        
        ze_diferente = False
        for b1, b2 in zip(mol_orig.GetBonds(), mol_sel.GetBonds()):
            if b1.GetBondType() == Chem.rdchem.BondType.DOUBLE and b2.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                if b1.GetStereo() != b2.GetStereo():
                    ze_diferente = True
                    break

        mol_mirror = Chem.Mol(mol_sel)
        for idx, conf in sel_dict.items():
            atom = mol_mirror.GetAtomWithIdx(idx)
            tag = atom.GetChiralTag()
            if tag == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
            elif tag == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
        Chem.AssignStereochemistry(mol_mirror, force=True, cleanIt=True)
        
        smi_mirror = Chem.MolToSmiles(mol_mirror, isomericSmiles=True)
        if smi_sel_iso == smi_mirror and not ze_diferente:
            return "identical"
            
        if ze_diferente:
            return "diastereomer"
            
        todos_invertidos = all(orig_dict.get(idx) != sel_dict.get(idx) for idx in orig_dict) if len(orig_dict) == len(sel_dict) and len(orig_dict) > 0 else False
        if todos_invertidos:
            return "enantiomer"
        else:
            return "diastereomer"
    except Exception:
        return "different"

def generar_multichoice_xml(questions_list, lang):
    """Genera el XML para Moodle de tipo 'multichoice' con las 4 opciones de respuesta exigidas."""
    t = TEXTS_STEREO[lang]
    quiz = ET.Element("quiz")
    
    opciones_map = {
        "enantiomer": t["rel_enantiomer"],
        "diastereomer": t["rel_diastereomer"],
        "identical": t["rel_identical"],
        "different": t["rel_different"]
    }
    
    for item in questions_list:
        name1 = item["name1"]
        name2 = item["name2"]
        correct_rel = item["relationship"]
        b64_img = item["image_b64"]
        
        question = ET.Element("question", type="multichoice")
        
        q_name = ET.SubElement(question, "name")
        ET.SubElement(ET.SubElement(q_name, "text"), "text").text = t["question_title"].format(name1, name2)
        
        q_text = ET.SubElement(question, "questiontext", format="html")
        img_tag = f'<p><img src="data:image/png;base64,{b64_img}" alt="Comparación estéreo" /></p>' if b64_img else ""
        ET.SubElement(q_text, "text").text = f"<![CDATA[<p>{t['question_text']}</p>{img_tag}]]>"
        
        for key, text_val in opciones_map.items():
            fraction = "100" if key == correct_rel else "0"
            answer = ET.SubElement(question, "answer", fraction=fraction, format="moodle_auto_format")
            ET.SubElement(answer, "text").text = text_val
            feedback = ET.SubElement(answer, "feedback", format="html")
            ET.SubElement(feedback, "text").text = "¡Correcto!" if fraction == "100" else "Incorrecto."
            
        quiz.append(question)
        
    tree = ET.ElementTree(quiz)
    ET.indent(tree, space="  ")
    return ET.tostring(tree.getroot(), encoding="utf-8", xml_declaration=True)

def render_nueva_herramienta(texts_main):
    """Función principal que renderiza el módulo de estereoisómeros en Streamlit."""
    lang = st.session_state.get("lang", "en")
    t = TEXTS_STEREO[lang]
    
    st.title(t["title"])
    st.write(t["desc"])
    
    if not RDKIT_AVAILABLE:
        st.warning(t["error_rdkit"])
        return

    # Inicializar estado para este módulo
    if "stereo_questions" not in st.session_state:
        st.session_state.stereo_questions = []

    # --- DISEÑO DE DOS COLUMNAS ---
    col_input, col_list = st.columns([1.1, 0.9], gap="large")

    # ================= COLUMNA IZQUIERDA: ENTRADA (INDIVIDUAL Y MASIVA) =================
    with col_input:
        tab_ind, tab_bulk = st.tabs([t["tab_individual"], t["tab_bulk"]])

        with tab_ind:
            # Paso 1: Introducir molécula base
            with st.form("form_base_stereo"):
                c_name = st.text_input(t["compound_name"])
                submitted_base = st.form_submit_button(t["btn_analyze"], type="primary")
                
                if submitted_base:
                    if not c_name.strip():
                        st.error("Por favor, introduce el nombre de una molécula." if lang == "es" else "Please enter a molecule name.")
                        st.session_state.active_mol_data = None
                    else:
                        with st.spinner("Buscando estructura y analizando elementos estereogénicos..."):
                            smiles = resolver_smiles(c_name)
                            if not smiles:
                                st.error(t["error_smiles"].format(c_name))
                                st.session_state.active_mol_data = None
                            else:
                                mol = Chem.MolFromSmiles(smiles)
                                Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
                                
                                # Comprobar centros quirales (incluyendo no asignados '?')
                                chiral_centers = FindMolChiralCenters(mol, includeUnassigned=True)
                                tiene_quirales_sin_asignar = any(config == '?' for idx, config in chiral_centers)
                                
                                # Comprobar dobles enlaces con estereoquímica Z/E explícita
                                has_double_bonds_specified = any(
                                    b.GetBondType() == Chem.rdchem.BondType.DOUBLE and 
                                    b.GetStereo() in [Chem.rdchem.BondStereo.STEREOZ, Chem.rdchem.BondStereo.STEREOE]
                                    for b in mol.GetBonds()
                                )
                                
                                if tiene_quirales_sin_asignar:
                                    st.warning("La molécula tiene centros quirales pero **no se ha indicado su estereoquímica** (ej: utiliza prefijos como *(2R)* o *(2S)*)." if lang == "es" else "The molecule has chiral centers but stereochemistry is not specified.")
                                    st.session_state.active_mol_data = None
                                elif not chiral_centers and not has_double_bonds_specified:
                                    st.warning("La molécula analizada no presenta elementos estereogénicos detectables o especificados." if lang == "es" else "The analyzed molecule has no detectable or specified stereogenic elements.")
                                    st.session_state.active_mol_data = None
                                else:
                                    mol_para_enum = Chem.Mol(mol)
                                    Chem.RemoveStereochemistry(mol_para_enum)
                                    options = StereoEnumerationOptions(unique=True)
                                    isomeros = list(EnumerateStereoisomers(mol_para_enum, options=options))
                                    
                                    lista_isomeros = []
                                    for idx_iso, isom in enumerate(isomeros, 1):
                                        desc = obtener_descripcion_estereo(isom)
                                        lista_isomeros.append({
                                            "id": idx_iso,
                                            "mol": isom,
                                            "config": desc
                                        })
                                        
                                    st.session_state.active_mol_data = {
                                        "name": c_name.strip(),
                                        "mol_orig": mol,
                                        "smiles_orig": smiles,
                                        "isomers": lista_isomeros
                                    }

            # Paso 2: Si la molécula base está cargada, permitir elegir modo de comparación
            if "active_mol_data" in st.session_state and st.session_state.active_mol_data:
                data = st.session_state.active_mol_data
                st.success(f"Molécula base cargada: `{data['name']}`")

                modo = st.radio(t["mode_selection"], [t["option_auto"], t["option_manual"]])

                textos_rel = {
                    "enantiomer": t["rel_enantiomer"],
                    "diastereomer": t["rel_diastereomer"],
                    "identical": t["rel_identical"],
                    "different": t["rel_different"]
                }

                if modo == t["option_auto"]:
                    if not data["isomers"]:
                        st.warning("La molécula base no generó estereoisómeros adicionales.")
                    else:
                        options_map_sel = {f"Estereoisómero {iso['id']} [{iso['config']}]": iso for iso in data["isomers"]}
                        seleccion_str = st.selectbox(t["select_isomer"], list(options_map_sel.keys()))
                        
                        if seleccion_str:
                            iso_elegido = options_map_sel[seleccion_str]
                            mol_sel = iso_elegido["mol"]
                            
                            relacion = calcular_relacion_estereoquimica(data["mol_orig"], mol_sel)
                            st.info(f"**{t['relationship_label']}** {textos_rel[relacion]}")
                            
                            img_b64 = generar_imagen_comparativa_b64(data["mol_orig"], mol_sel)
                            if img_b64:
                                st.image(f"data:image/png;base64,{img_b64}", width=350)
                                
                            if st.button(t["add_question"], type="primary", icon=":material/add_task:", use_container_width=True):
                                st.session_state.stereo_questions.append({
                                    "name1": f"{data['name']} (Orig)",
                                    "name2": f"Isómero {iso_elegido['id']} [{iso_elegido['config']}]",
                                    "relationship": relacion,
                                    "image_b64": img_b64
                                })
                                st.success(t["success_added"])

                else:  # Modo manual para el segundo compuesto
                    with st.form("form_manual_comparison"):
                        comp2_input = st.text_input(t["manual_comparison_name"])
                        submitted_eval = st.form_submit_button(t["btn_evaluate_manual"], type="primary")
                        
                        if submitted_eval:
                            if not comp2_input.strip():
                                st.error("Por favor, introduce un nombre o SMILES válido.")
                                st.session_state.manual_eval_result = None
                            else:
                                smi2 = resolver_smiles(comp2_input)
                                if not smi2:
                                    st.error(t["error_smiles"].format(comp2_input))
                                    st.session_state.manual_eval_result = None
                                else:
                                    mol2 = Chem.MolFromSmiles(smi2)
                                    relacion = calcular_relacion_estereoquimica(data["mol_orig"], mol2)
                                    img_b64 = generar_imagen_comparativa_b64(data["mol_orig"], mol2)
                                    
                                    st.session_state.manual_eval_result = {
                                        "name2": comp2_input.strip(),
                                        "relationship": relacion,
                                        "image_b64": img_b64
                                    }

                    # Mostrar resultado de la evaluación manual fuera del form
                    if "manual_eval_result" in st.session_state and st.session_state.manual_eval_result:
                        res = st.session_state.manual_eval_result
                        st.info(f"**{t['relationship_label']}** {textos_rel[res['relationship']]}")
                        if res["image_b64"]:
                            st.image(f"data:image/png;base64,{res['image_b64']}", width=350)
                            
                        if st.button(t["add_question"], type="primary", icon=":material/add_task:", use_container_width=True, key="add_manual_q"):
                            st.session_state.stereo_questions.append({
                                "name1": f"{data['name']} (Orig)",
                                "name2": res["name2"],
                                "relationship": res["relationship"],
                                "image_b64": res["image_b64"]
                            })
                            st.success(t["success_added"])

        with tab_bulk:
            st.info(t["bulk_note"])
            uploaded_file = st.file_uploader(t["bulk_uploader"], type=["xlsx", "csv"], key="bulk_stereo_file")
            
            if uploaded_file and st.button(t["start_bulk"]):
                try:
                    if uploaded_file.name.endswith(".csv"):
                        df = pd.read_csv(uploaded_file)
                    else:
                        df = pd.read_excel(uploaded_file)
                        
                    if not {"molecule1", "molecule2"}.issubset(df.columns):
                        st.error("El archivo debe contener las columnas exactas: 'molecule1' y 'molecule2'.")
                    else:
                        success_count = 0
                        fail_count = 0
                        progress_bar = st.progress(0, text="Procesando lote de pares...")
                        total = len(df)
                        
                        for idx, row in df.iterrows():
                            m1_input = str(row["molecule1"]).strip()
                            m2_input = str(row["molecule2"]).strip()
                            progress_bar.progress((idx + 1) / total, text=f"Procesando: {m1_input} vs {m2_input}")
                            
                            smi1 = resolver_smiles(m1_input)
                            smi2 = resolver_smiles(m2_input)
                            
                            if smi1 and smi2:
                                mol1 = Chem.MolFromSmiles(smi1)
                                mol2 = Chem.MolFromSmiles(smi2)
                                relacion = calcular_relacion_estereoquimica(mol1, mol2)
                                img_b64 = generar_imagen_comparativa_b64(mol1, mol2)
                                
                                st.session_state.stereo_questions.append({
                                    "name1": m1_input,
                                    "name2": m2_input,
                                    "relationship": relacion,
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
        
        if st.session_state.stereo_questions:
            st.write(f"Total de preguntas en la lista: **{len(st.session_state.stereo_questions)}**")
            
            # --- BOTONES PRINCIPALES EN VERTICAL ---
            xml_data = generar_multichoice_xml(st.session_state.stereo_questions, lang)
            
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
                file_name="moodle_stereoisomers.xml",
                icon=":material/download:",
                mime="application/xml",
                use_container_width=True
            )

            if st.button(t["clear_all"], type="secondary", icon=":material/delete_sweep:", use_container_width=True):
                st.session_state.stereo_questions = []
                st.rerun()
            
            st.divider()
                
            textos_rel_list = {
                "enantiomer": t["rel_enantiomer"],
                "diastereomer": t["rel_diastereomer"],
                "identical": t["rel_identical"],
                "different": t["rel_different"]
            }
            
            # Contenedor con scroll para la lista de preguntas
            container_list = st.container(height=420)
            with container_list:
                for i, q in enumerate(st.session_state.stereo_questions):
                    cols_item = st.columns([1, 2.5, 0.6])
                    with cols_item[0]:
                        if q['image_b64']:
                            st.image(f"data:image/png;base64,{q['image_b64']}", width=100)
                    with cols_item[1]:
                        rel_texto = textos_rel_list.get(q['relationship'], q['relationship'])
                        st.markdown(f"**#{i+1}**<br>`{q['name1']}`<br>vs `{q['name2']}`<br>Tipo: **{rel_texto}**", unsafe_allow_html=True)
                    with cols_item[2]:
                        if st.button("🗑️", key=f"del_stereo_{i}"):
                            st.session_state.stereo_questions.pop(i)
                            st.rerun()
                    st.divider()
        else:
            st.info(t["empty_list"])