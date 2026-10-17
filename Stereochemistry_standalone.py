# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 09:05:18 2026

@author: Carlos Fernandez
"""

"""
Script de prueba independiente para detectar y resaltar centros quirales con RDKit.
"""

import sys
try:
    from rdkit import Chem
    from rdkit.Chem import FindMolChiralCenters, Draw
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
    Busca la molécula ingresada, calcula su estereoquímica, identifica múltiples centros quirales,
    los numera explícitamente y genera la visualización interactiva.
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
    
    # 4. Encontrar centros quirales (devuelve lista de tuplas: (índice_átomo, configuración))
    centros_quirales = FindMolChiralCenters(mol, includeUnassigned=True)
    
    indices_quirales = []
    notacion_centros = []
    
    if not centros_quirales:
        print("Aviso: La molécula cargada no posee centros quirales detectados o no están definidos espacialmente.")
    else:
        print(f"Centros quirales detectados (Índice, Configuración): {centros_quirales}")
        for idx, config in centros_quirales:
            indices_quirales.append(idx)
            notacion_centros.append(f"C{idx}:{config}")
            print(f" -> Carbono quiral en índice [{idx}] con configuración absoluta: ({config})")
            
            # Asignar una nota personalizada al átomo para que RDKit la dibuje junto al carbono
            # Ejemplo: C1(S), C2(R)
            atom = mol.GetAtomWithIdx(idx)
            atom.SetProp('atomNote', f"{idx}")
        
        formato_moodle_multiple = ", ".join(notacion_centros)
        print(f"\n💡 [Propuesta Moodle] Formato combinado para múltiples centros: # {formato_moodle_multiple}")

    # 5. Configurar opciones de dibujo de RDKit
    opts = Draw.MolDrawOptions()
    opts.addAtomIndices = True  # Muestra los índices de todos los átomos
    opts.fixedFontSize = 14     # Tamaño de etiqueta legible
    
    # 6. Generar imagen PIL de la molécula
    print("Renderizando imagen molecular...")
    try:
        img = Draw.MolToImage(
            mol, 
            size=(700, 550), 
            highlightAtoms=indices_quirales,
            drawOptions=opts
        )
    except Exception as e:
        print(f"Error al generar la imagen gráfica de la molécula: {e}")
        return
    
    # 7. Mostrar la imagen usando matplotlib
    fig = None
    try:
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.imshow(img)
        ax.axis('off')
              
        print("Abriendo ventana gráfica... (Cierra la ventana para continuar)")
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
    print("=== Generador de Estereoquímica RDKit ===")
    while True:
        try:
            nombre_usuario = input("\nIntroduce el nombre de una molécula en inglés (o 'salir'): ")
            if nombre_usuario.strip().lower() in ['salir', 'exit', 'q']:
                print("¡Hasta luego!")
                break
            if not nombre_usuario.strip():
                continue
            analizar_y_mostrar_quiralidad(nombre_usuario)
        except KeyboardInterrupt:
            print("\n¡Hasta luego!")
            break