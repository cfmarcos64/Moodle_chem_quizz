# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 09:05:18 2026

@author: Carlos Fernandez
"""

"""
Script independiente para detectar centros quirales, enumerar
estereoisómeros, comparar visualmente y clasificar la relación estereoquímica.
"""

import sys
try:
    from rdkit import Chem
    from rdkit.Chem import FindMolChiralCenters, Draw
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
    import requests
    import matplotlib
    import matplotlib.pyplot as plt
    
    # Configurar backend de matplotlib si es necesario
    try:
        if 'matplotlib.backends' not in sys.modules:
            matplotlib.use('QtAgg', force=False)
    except Exception:
        pass
        
    print("Todas las librerías necesarias se importaron correctamente.")
except ImportError as e:
    print(f"Error: Falta alguna librería requerida ({e}). Asegúrate de tener RDKit, requests, matplotlib y pillow.")
    sys.exit(1)

def buscar_smiles_nci(nombre_molecula):
    """
    Busca el SMILES canónico de una molécula a partir de sus nombres comunes o IUPAC
    utilizando el servicio web de NCI/CADD Chemical Identifier Resolver (CIR) con timeout.
    """
    nombre_limpio = nombre_molecula.strip()
    print(f"Buscando '{nombre_limpio}' en la API de NCI CIR (esto puede tardar unos segundos)...")
    
    url = f"http://cactus.nci.nih.gov/chemical/structure/{requests.utils.quote(nombre_limpio)}/smiles"
    try:
        respuesta = requests.get(url, timeout=8)
        if respuesta.status_code == 200:
            smiles = respuesta.text.strip().split('\n')[0]
            if smiles and "Error" not in smiles:
                print(f" -> SMILES encontrado: {smiles}")
                return smiles
    except requests.exceptions.Timeout:
        print("Error: La solicitud a la API de NCI CIR ha tardado demasiado (Timeout).")
    except Exception as ex:
        print(f"Error de conexión con la API: {ex}.")
        
    print(f"No se pudo encontrar un SMILES válido para '{nombre_molecula}'.")
    return None

def analizar_y_mostrar_quiralidad(nombre_molecula):
    """
    Busca la molécula ingresada, calcula su estereoquímica, identifica centros quirales,
    enumera los estereoisómeros, clasifica la relación estereoquímica y genera la comparativa visual limpia.
    """
    print("\n==========================================")
    print(f" Analizando compuesto: {nombre_molecula}")
    print("==========================================")
    
    plt.close('all')
    
    # 1. Obtener SMILES desde la API de NCI CIR
    smiles = buscar_smiles_nci(nombre_molecula)
    if not smiles:
        return
        
    # 2. Cargar la molécula en RDKit
    print("Generando estructura molecular...")
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("Error: RDKit no pudo parsear el SMILES obtenido.")
        return
        
    # 3. Forzar cálculo y asignación de estereoquímica (etiquetas CIP R/S)
    try:
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    except Exception as e:
        print(f"Aviso al asignar estereoquímica: {e}")
    
    # 4. Encontrar centros quirales del original
    centros_quirales = FindMolChiralCenters(mol, includeUnassigned=True)
    
    if not centros_quirales:
        print("Aviso: La molécula cargada no posee centros quirales detectados o no están definidos espacialmente.")
        return
    else:
        print(f"Centros quirales detectados (Índice, Configuración): {centros_quirales}")
        notacion_centros = [f"C{idx}:{config}" for idx, config in centros_quirales]
        formato_moodle_multiple = ", ".join(notacion_centros)
        print(f"\n💡 [Propuesta Moodle] Formato combinado para múltiples centros: # {formato_moodle_multiple}")

        # 5. Enumerar todas las combinaciones posibles de estereoisómeros
        print("\n------------------------------------------")
        print(" 🔄 Enumerando todas las posibilidades de estereoisómeros:")
        print("------------------------------------------")
        try:
            mol_para_enum = Chem.Mol(mol)
            Chem.RemoveStereochemistry(mol_para_enum)
            
            options = StereoEnumerationOptions(unique=True)
            isomeros = list(EnumerateStereoisomers(mol_para_enum, options=options))
            print(f"Se han encontrado {len(isomeros)} estereoisómero(s) único(s) para esta estructura:")
            
            lista_isomeros_obj = []
            for i, isom in enumerate(isomeros, 1):
                Chem.AssignStereochemistry(isom, force=True, cleanIt=True)
                centros_isom = FindMolChiralCenters(isom, includeUnassigned=True)
                smi_isom = Chem.MolToSmiles(isom)
                config_str = ", ".join([f"C{idx}:{conf}" for idx, conf in centros_isom])
                print(f"  [{i}] Estereoisómero {i}: [{config_str}] -> SMILES: {smi_isom}")
                lista_isomeros_obj.append((isom, config_str))
        except Exception as e:
            print(f"No se pudieron enumerar los estereoisómeros automáticamente: {e}")
            return

    # 6. Pedir al usuario que elija un estereoisómero de la lista
    while True:
        try:
            seleccion = input(f"\nElige el número de estereoisómero a comparar (1 a {len(lista_isomeros_obj)}): ")
            idx_elegido = int(seleccion) - 1
            if 0 <= idx_elegido < len(lista_isomeros_obj):
                mol_seleccionado, config_elegida = lista_isomeros_obj[idx_elegido]
                break
            else:
                print("Número fuera de rango. Inténtalo de nuevo.")
        except ValueError:
            print("Por favor, introduce un número válido.")

    # 7. Calcular y clasificar la relación estereoquímica
    orig_dict = dict(centros_quirales)
    sel_dict = dict(FindMolChiralCenters(mol_seleccionado, includeUnassigned=True))
    
    # Comprobar si el seleccionado es una forma MESO
    mol_mirror = Chem.Mol(mol_seleccionado)
    for idx, conf in sel_dict.items():
        atom = mol_mirror.GetAtomWithIdx(idx)
        tag = atom.GetChiralTag()
        if tag == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
            atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif tag == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
            atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
    Chem.AssignStereochemistry(mol_mirror, force=True, cleanIt=True)
    
    smi_sel = Chem.MolToSmiles(mol_seleccionado, isomericSmiles=True)
    smi_mirror = Chem.MolToSmiles(mol_mirror, isomericSmiles=True)
    es_meso = (smi_sel == smi_mirror)
    
    # Determinar el tipo de relación
    print("\n------------------------------------------")
    print(" 🔬 Análisis Estereoquímico:")
    print("------------------------------------------")
    if orig_dict == sel_dict:
        relacion = "Es la misma molécula (idéntica al original)."
    elif es_meso:
        relacion = "Es una forma MESO (molécula aquiral globalmente a pesar de poseer centros quirales)."
    else:
        todos_invertidos = all(orig_dict[idx] != sel_dict.get(idx) for idx in orig_dict)
        if todos_invertidos:
            relacion = "Son ENANTIÓMEROS (imágenes especulares no superponibles)."
        else:
            relacion = "Son DIASTEREOISÓMEROS (estereoisómeros que no son imágenes especulares)."
            
    print(f" -> Relación detectada: {relacion}")

    # 8. Configurar opciones de dibujo limpias (sin índices de átomos)
    opts = Draw.MolDrawOptions()
    opts.addAtomIndices = False

    # 9. Generar imagen comparativa limpia en cuadrícula (Original vs Seleccionado)
    print("\nRenderizando imagen comparativa limpia...")
    try:
        img = Draw.MolsToGridImage(
            [mol, mol_seleccionado],
            molsPerRow=2,
            subImgSize=(450, 450),
            legends=None,
            drawOptions=opts
        )
    except Exception as e:
        print(f"Error al generar la imagen comparativa: {e}")
        return
    
    # 10. Mostrar la imagen usando matplotlib
    fig = None
    try:
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.imshow(img)
        ax.axis('off')
              
        print("Abriendo ventana gráfica comparativa... (Cierra la ventana para continuar)")
        plt.show(block=False)
        
        while plt.fignum_exists(fig.number):
            plt.pause(0.2)
            
        print("Ventana cerrada por el usuario. Continuando...")
    except KeyboardInterrupt:
        print("\nInterrupción detectada. Cerrando ventana gráfica...")
        if fig is not None:
            plt.close(fig)
    except Exception as e:
        print(f"Error al mostrar la ventana gráfica: {e}")
    finally:
        if fig is not None:
            plt.close(fig)
        plt.close('all')

if __name__ == "__main__":
    print("=== Generador y Clasificador de Estereoisómeros RDKit ===")
    while True:
        try:
            nombre_usuario = input("\nIntroduce el nombre de una molécula en inglés (o 'salir'): ")
            if nombre_usuario.strip().lower() in ['salir', 'exit', 'q']:
                print("¡Hasta luegokum!" if '¡Hasta luego!' else "¡Hasta luego!")
                break
            if not nombre_usuario.strip():
                continue
            analizar_y_mostrar_quiralidad(nombre_usuario)
        except KeyboardInterrupt:
            print("\n¡Hasta luego!")
            break