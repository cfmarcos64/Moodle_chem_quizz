# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 09:05:18 2026

@author: Carlos F. Marcos
"""

import streamlit as st

# --- 1. PAGE CONFIGURATION ---
st.set_page_config(
    page_title="Moodle Chemistry Suite",
    layout="wide",
    page_icon="🧪",
    initial_sidebar_state="expanded" 
)

# --- 2. MODULE IMPORTATION ---
try:
    from Nomenclature_for_Moodle import render_nueva_herramienta as render_nomenclature_app
    from FG_for_Moodle import render_nueva_herramienta as render_fg_app
    from Stereochemistry_for_Moodle import render_nueva_herramienta as render_stereo_app
    from Compare_for_Moodle import render_nueva_herramienta as render_compare_app
except ImportError as e:
    st.error(f"Error cargando módulos de la suite: {e}")

# --- 3. LANGUAGE CENTRAL MANAGEMENT ---
if "lang" not in st.session_state:
    st.session_state.lang = "en"

def toggle_language():
    if st.session_state.lang_toggle:
        st.session_state.lang = "es"
    else:
        st.session_state.lang = "en"
    
# Multilingual texts for the navigation suite
ui_texts = {
    "es": {
        "sidebar_title": "🧪 Suite de Química para Moodle", 
        "selector_label": "Selecciona una herramienta:", 
        "opt": [
            "Nomenclatura Sistemática", 
            "Grupos Funcionales", 
            "Determinación de Estereoquímica", 
            "Comparación de Moléculas"
        ]
    },
    "en": {
        "sidebar_title": "🧪 Moodle Chemistry Suite", 
        "selector_label": "Select a tool:", 
        "opt": [
            "Systematic Nomenclature", 
            "Functional Groups", 
            "Stereochemistry", 
            "Molecule Comparison"
        ]
    }
}

# --- 4. SIDEBAR ---

with st.sidebar:
    # Language toggle UI
    style_es = "color: #28a745; font-weight: bold;" if st.session_state.lang == "es" else "color: #888;"
    style_en = "color: #28a745; font-weight: bold;" if st.session_state.lang == "en" else "color: #888;"
    
    col_txt, col_tgl = st.columns([0.7, 0.3])
    
    with col_txt:
        st.markdown(
            f'''
            <div style="
                background-color: rgba(255, 255, 255, 0.4); 
                padding: 5px 10px; 
                border-radius: 15px; 
                text-align: center;
                margin-top: 5px;
                font-size: 0.85rem;
                display: flex;
                justify-content: center;
                gap: 10px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.05);
            ">
                <span style="{style_en}">ENGLISH</span>
                <span style="color: #ccc;">|</span>
                <span style="{style_es}">ESPAÑOL</span>
            </div>
            ''', 
            unsafe_allow_html=True
        )

    with col_tgl:
        st.toggle(
            " ", 
            value=(st.session_state.lang == "es"), 
            key="lang_toggle", 
            on_change=toggle_language,
            label_visibility="collapsed"
        )
    
    st.divider()
    
    # Title
    st.title(ui_texts[st.session_state.lang]["sidebar_title"])
    
    # Tool selector
    opcion = st.radio(
        ui_texts[st.session_state.lang]["selector_label"], 
        ui_texts[st.session_state.lang]["opt"]
    )

# --- Colors dictionary for tools ---
tool_colors = {
    "Nomenclatura": "#e3f2fd",   # Light blue
    "Nomenclature": "#e3f2fd",   
    "Grupos": "#f1f8e9",         # Light green
    "Functional": "#f1f8e9",     
    "Determinación": "#fff3e0",# Light orange
    "Stereochemistry": "#fff3e0",  
    "Comparación": "#fce4ec",    # Light pink
    "Comparison": "#fce4ec"      
}   

# --- 5. ROUTING ---

def set_sidebar_color(option):
    color = "#ffffff" # Default white
    for key, hex_color in tool_colors.items():
        if key in option:
            color = hex_color
            break
            
    st.markdown(f"""
        <style>
            [data-testid="stSidebar"] {{
                background-color: {color};
            }}
        </style>
    """, unsafe_allow_html=True)

set_sidebar_color(opcion)

# Language-independent option routing
opciones_en = ui_texts["en"]["opt"]
opciones_es = ui_texts["es"]["opt"]

if opcion == opciones_en[0] or opcion == opciones_es[0]:
    render_nomenclature_app(ui_texts)
elif opcion == opciones_en[1] or opcion == opciones_es[1]:
    render_fg_app(ui_texts)
elif opcion == opciones_en[2] or opcion == opciones_es[2]:
    render_stereo_app(ui_texts)
elif opcion == opciones_en[3] or opcion == opciones_es[3]:
    render_compare_app(ui_texts)