# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 17:14:25 2024

@author: giuli
"""


import os
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import subprocess
from glob import glob

def check_tool(tool_name):
    """Verifica si un programa está instalado y disponible en la terminal."""
    try:
        subprocess.run([tool_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except FileNotFoundError:
        return False
    except subprocess.CalledProcessError:
        return True  # El programa existe pero no acepta '--version'

def read_csv_files():
    """Lee todos los archivos CSV en el mismo directorio donde se encuentra el script y los almacena en un diccionario."""
    script_dir = os.path.dirname(os.path.abspath(__file__))  # Obtiene la ruta del script
    files = glob(os.path.join(script_dir, "*.csv"))
    
    if not files:
        print("❌ No se encontraron archivos CSV en:", script_dir)
    
    print(f"📂 Archivos encontrados: {files}")  # Debugging
    
    dataframes = {}
    consensus_sequences = {}
    for file in files:
        df = pd.read_csv(file)
        
        # Asegurar que la columna 'Posición' es el índice
        if 'Posición' in df.columns:
            df.set_index('Posición', inplace=True)
        
        # Extraer la secuencia consenso
        if 'Residuo Consenso' in df.columns:
            consensus_sequences[os.path.basename(file)] = "".join(df['Residuo Consenso'].tolist())
            df = df.drop(columns=['Residuo Consenso'])  # Eliminar esta columna antes de procesar
        
        # Convertir a valores numéricos
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
        
        # Calcular la información total por posición y asegurar que logomaker no normalice
        if 'Información (bits)' in df.columns:
            total_info = df['Información (bits)']
            df = df.drop(columns=['Información (bits)'])
            df = df.multiply(total_info, axis=0)
        
        dataframes[os.path.basename(file)] = df
    
    return dataframes, consensus_sequences

def align_consensus_sequences(consensus_sequences):
    """Alinea las secuencias consenso usando MUSCLE o ClustalO y devuelve la alineación."""
    alignment_file = "consensus_alignment.fasta"
    aligned_file = "consensus_aligned.fasta"
    
    with open(alignment_file, "w") as f:
        for name, seq in consensus_sequences.items():
            f.write(f">{name}\n{seq}\n")
    
    # Determinar qué herramienta usar
    if check_tool("muscle"):
        cmd = ["muscle", "-in", alignment_file, "-out", aligned_file]
    elif check_tool("clustalo"):
        cmd = ["clustalo", "-i", alignment_file, "-o", aligned_file, "--force"]
    else:
        raise RuntimeError("❌ No se encontró MUSCLE ni Clustal Omega. Instálalos antes de ejecutar el script.")
    
    subprocess.run(cmd, check=True)
    
    aligned_sequences = {}
    aligned_order = []
    with open(aligned_file, "r") as f:
        name = ""
        seq = ""
        for line in f:
            if line.startswith(">"):
                if name:
                    aligned_sequences[name] = seq
                    aligned_order.append(name)
                name = line.strip()[1:]
                seq = ""
            else:
                seq += line.strip()
        if name:
            aligned_sequences[name] = seq
            aligned_order.append(name)
    
    return aligned_sequences, aligned_order

def reorder_dataframe(df, aligned_seq):
    """Reordena el DataFrame de frecuencias según la alineación de la secuencia consenso."""
    new_index = list(range(len(aligned_seq)))
    new_df = pd.DataFrame(0, index=new_index, columns=df.columns)
    
    original_positions = df.index.tolist()
    pos_map = {i: original_positions[i] for i in range(len(original_positions))}
    
    mapped_positions = []
    current_pos = 0
    for i, aa in enumerate(aligned_seq):
        if aa != "-" and current_pos < len(original_positions):
            mapped_positions.append(pos_map[current_pos])
            current_pos += 1
        else:
            mapped_positions.append(None)
    
    new_df.index = new_index
    
    for i, pos in enumerate(mapped_positions):
        if pos is not None and pos in df.index:
            new_df.loc[i] = df.loc[pos]
    
    return new_df.fillna(0)

def align_sequence_logos(dataframes, aligned_sequences, aligned_order, output_file):
    """Genera logos alineados según la alineación de las secuencias consenso y los exporta como PNG y SVG."""
    fig, axes = plt.subplots(len(dataframes), 1, figsize=(10, 3 * len(dataframes)), sharex=True)
    
    if len(dataframes) == 1:
        axes = [axes]  # Convertir en lista si solo hay un logo
    
    for ax, name in zip(axes, aligned_order):
        df = reorder_dataframe(dataframes[name], aligned_sequences[name])
        logomaker.Logo(df, ax=ax, font_name='Arial Rounded MT Bold', color_scheme='chemistry')
        ax.set_xlim(-0.5, len(df)-0.5)
        ax.set_title(name)
        ax.set_ylabel("Información (bits)")
        ax.set_ylim(0, None)  # Ajusta el eje Y dinámicamente
    
    plt.xlabel("Posición alineada en la secuencia")
    plt.savefig(output_file, dpi=300)
    plt.savefig(output_file.replace(".png", ".svg"), format="svg")  # Guardar como SVG
    
    plt.show()

# Uso del programa
dataframes, consensus_sequences = read_csv_files()
aligned_sequences, aligned_order = align_consensus_sequences(consensus_sequences)
align_sequence_logos(dataframes, aligned_sequences, aligned_order, "aligned_logos.png")



