# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 17:14:25 2024

@author: giuli
"""


from Bio import AlignIO
import math
import pandas as pd

# Cargar el alineamiento en formato FASTA
alignment = AlignIO.read("Seq_Consenso_Clusters_Johnma_FASTA_MSA_exported_all_UPGAMA_extremosRecortados.fas", "fasta")

# Crear listas para almacenar frecuencias, entropía y la secuencia consenso
frecuencias = []
informacion = []
consenso = []

# Entropía máxima considerando solo los 20 aminoácidos estándar
H_max = math.log2(20)

# Iterar sobre cada posición en el alineamiento
for i in range(alignment.get_alignment_length()):
    columna = [seq[i] for seq in alignment]
    
    # Filtrar los gaps ('-') de la columna
    columna_sin_gaps = [res for res in columna if res != '-']
    total_residuos = len(columna_sin_gaps)
    
    # Evitar cálculos si la columna está completamente vacía (solo gaps)
    if total_residuos == 0:
        frecuencias.append({})
        informacion.append(None)  # Información no definida para posiciones sin datos
        consenso.append('-')  # Gap para posiciones vacías
        continue
    
    # Calcular las frecuencias de cada residuo
    conteo = {res: columna_sin_gaps.count(res) / total_residuos for res in set(columna_sin_gaps)}
    frecuencias.append(conteo)
    
    # Calcular la entropía de Shannon
    entropia = -sum(freq * math.log2(freq) for freq in conteo.values())
    
    # Calcular la información
    informacion.append(H_max - entropia)

    # Determinar el residuo consenso (máxima frecuencia)
    residuo_consenso = max(conteo, key=conteo.get)
    consenso.append(residuo_consenso)

# Crear un DataFrame para guardar los resultados
df_frecuencias = pd.DataFrame(frecuencias).fillna(0)  # Completa valores ausentes con 0
df_frecuencias["Información (bits)"] = informacion
df_frecuencias["Residuo Consenso"] = consenso

# Guardar el resultado como CSV
df_frecuencias.to_csv("frecuencias_informacion_sin_gaps.csv", index_label="Posición")

# Generar la secuencia consenso completa
secuencia_consenso = ''.join(consenso)

# Guardar la secuencia consenso en un archivo separado
with open("secuencia_consenso.txt", "w") as f:
    f.write(f">Secuencia_Consenso\n{secuencia_consenso}\n")

print("Tabla exportada a frecuencias_informacion_sin_gaps.csv")
print("Secuencia consenso exportada a secuencia_consenso.txt")
