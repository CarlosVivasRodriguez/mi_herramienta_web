import os
import zipfile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import requests
from goatools.obo_parser import GODag
from pathlib import Path
import shutil
import re
import numpy as np
from flask import Flask, render_template, request, redirect, url_for, flash, send_file, session, jsonify  # Mueve jsonify aquí
from matplotlib import colormaps
from matplotlib.gridspec import GridSpec  # Para las subgráficas en figura 3
from matplotlib.ticker import MaxNLocator  # Asegúrate de importar esto
from itertools import cycle  # Para generar combinaciones de colores
from upsetplot import UpSet, from_memberships  # Para generar UpSet plots
import warnings  # Para ignorar advertencias que no necesitas en el contexto de la app
import io  # Para manejo de buffers en memoria
import matplotlib.colors as mcolors
import webbrowser  # Importa el módulo para abrir el navegador automáticamente
import threading  # Para manejar el temporizador que abre el navegador con retraso
import matplotlib.cm as cm
import time
import gzip
from goatools.anno.gaf_reader import GafReader
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudy
import traceback 

# Importamos y configuramos matplotlib para evitar errores de GUI al generar figuras
import matplotlib
matplotlib.use('Agg') 
# Inicializamos la aplicación Flask aquí
app = Flask(__name__)

# Ignorar advertencias
warnings.simplefilter(action='ignore', category=FutureWarning)

app.secret_key = 'your_secret_key'

# Inicializar las variables globales para su uso en toda la aplicación
orthogroups_df = None
gene_count_df = None
species_selected = None
ortogrupos_por_combinacion = None
abreviaturas = None
selected_orthogroups = None
species_urls = None 
selected_orthogroups_list = None
final_protein_list = None

def cargar_datos_carpeta(ruta_carpeta):
    print(f"Archivos en la carpeta '{ruta_carpeta}':")
    print(os.listdir(ruta_carpeta))
    
    gene_count_path = os.path.join(ruta_carpeta, 'Orthogroups.GeneCount.tsv')
    orthogroups_path = os.path.join(ruta_carpeta, 'Orthogroups.tsv')
    single_copy_path = os.path.join(ruta_carpeta, 'Orthogroups_SingleCopyOrthologues.txt')
    unassigned_path = os.path.join(ruta_carpeta, 'Orthogroups_UnassignedGenes.tsv')

    if not os.path.exists(gene_count_path):
        raise FileNotFoundError(f"No se encontró el archivo: {gene_count_path}")
    if not os.path.exists(orthogroups_path):
        raise FileNotFoundError(f"No se encontró el archivo: {orthogroups_path}")
    if not os.path.exists(single_copy_path):
        raise FileNotFoundError(f"No se encontró el archivo: {single_copy_path}")
    if not os.path.exists(unassigned_path):
        raise FileNotFoundError(f"No se encontró el archivo: {unassigned_path}")

    gene_count_df = pd.read_csv(gene_count_path, sep='\t')
    orthogroups_df = pd.read_csv(orthogroups_path, sep='\t')
    single_copy_df = pd.read_csv(single_copy_path, sep='\t')
    unassigned_df = pd.read_csv(unassigned_path, sep='\t')
    
    return gene_count_df, orthogroups_df, single_copy_df, unassigned_df

def eliminar_archivo(ruta_archivo):
    if os.path.exists(ruta_archivo):
        os.remove(ruta_archivo)

def crear_directorio_plots():
    if not os.path.exists(os.path.join('static', 'plots')):
        os.makedirs(os.path.join('static', 'plots'))

def generar_figura_1(gene_count_df):
    fig, ax = plt.subplots(figsize=(10, 6))

 # Obtener los datos y reemplazar valores infinitos o nulos
    data = gene_count_df['Total'].replace([np.inf, -np.inf], np.nan).dropna()

    # Crear bins personalizados (ajustando correctamente los límites para cada rango)
    bins = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 15.5, 18.5, 20.5, 25.5, 30.5, 35.5, 
        40.5, 50.5, 60.5, 70.5, 80.5, 90.5, 100.5, 200.5, 300.5, 400.5, 500.5, 1000.5, int(data.max())+1]

    # Crear etiquetas personalizadas para los grupos
    labels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11-12', '13-15', '16-18', '19-20', '21-25', 
          '26-30', '31-35', '36-40', '41-50', '51-60', '61-70', '71-80', '81-90', '91-100', '101-200', 
          '201-300', '301-400', '401-500', '501-1000', '1001+']


    # Agrupar los datos en los bins definidos
    grouped_data = pd.cut(data, bins=bins, labels=labels, include_lowest=True, right=False)

    # Contar las frecuencias en cada bin
    counts = grouped_data.value_counts().sort_index()  # CORREGIDO: esta línea ahora está correctamente indentada

    # Crear gráfico de barras con los datos agrupados
    sns.barplot(x=counts.index, y=counts.values, ax=ax, color='orange')

    # Añadir líneas verticales separadoras entre los grupos solicitados
    separator_positions = {
        '10_11-12': 9.5,   # Entre 10 y 11-12
        '19-20_21-25': 13.5,  # Entre 19-20 y 21-25
        '36-40_41-50': 17.5,  # Entre 36-40 y 41-50
        '91-100_101-200': 23.5,  # Entre 91-100 y 101-200
        '401-500_501-1000': 27.5  # Entre 401-500 y 501-1000
    }

    # Dibujar las líneas verticales en las posiciones especificadas
    for position in separator_positions.values():
        ax.axvline(x=position, color='black', linestyle='--')

    # Añadir explicaciones dentro del gráfico para cada rango
    text_positions = {
        (5, 5): 'INC by 1 ',      # Para el rango de 1 a 10
        (10, 13): 'INC by 2 ',   # Para el rango de 11-12
        (15.5, 15.5): 'INC by 5 ',   # Para el rango de 21-25
        (20.5, 20.5): 'INC by 10 ',  # Para el rango de 41-50
        (25.5, 25.5): 'INC by 100 ', # Para el rango de 101-200
        (28.5, 28.5): 'INC by 500 ', # Para el rango de 501-1000
    }
    # Altura fija para todos los textos (ligeramente por encima de la barra más alta del gráfico)
    fixed_height = counts.max() * 0.7  # Puedes ajustar este valor para mayor altura
    # Colocar los textos explicativos dentro del gráfico en las posiciones deseadas
    for (start, end), text in text_positions.items():
        # Calcular el midpoint (punto medio) entre las posiciones de las etiquetas
        midpoint = (start + end) / 2

        # Colocar el texto en vertical en una altura fija
        ax.text(midpoint, fixed_height, text, ha='center', fontsize=14, color='grey', rotation='vertical')

    # Títulos y etiquetas
    ax.set_title('Protein Distribution by Orthogroup', fontsize=23)
    ax.set_xlabel('Number of Proteins', fontsize=20)
    ax.set_ylabel('Frequency', fontsize=20)

    # Ajustar las etiquetas del eje X para que sean visibles
    ax.tick_params(axis='x', rotation=90)

    # Ajustar los límites del eje Y para acomodar las etiquetas encima de las barras
    ax.set_ylim(0, counts.max() + 220)

    # Añadir etiquetas encima de cada barra con el número de ortogrupos
    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.text(p.get_x() + p.get_width() / 2, height + 0.8, f'{int(height)}',
                    ha='center', va='bottom', fontsize=12)

    # Ajustar el espaciado para que el gráfico llene mejor el área disponible
    plt.tight_layout()

    crear_directorio_plots()
    
    image_path = os.path.join('static', 'plots', 'figura_1.png')
    plt.savefig(image_path)
    plt.close(fig)

    return 'figura_1.png'

def generar_figura_2(gene_count_df, axes=None):
    # Crear los ejes si no se proporcionan
    if axes is None:
        fig, axes = plt.subplots(figsize=(10, 6))

    # Convertir todas las columnas numéricas excepto la primera (ID de ortogrupos)
    gene_count_df.iloc[:, 1:] = gene_count_df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')

    # Seleccionar las especies (todas las columnas excepto la primera y la última)
    species = gene_count_df.columns[1:-1]
    
    # Contar cuántas especies comparten ortogrupos
    shared_orthogroups = gene_count_df[species].apply(lambda x: x > 0).sum(axis=1)
    unique_orthogroups = shared_orthogroups.value_counts()

    # Luego generar el gráfico como antes
    sns.barplot(x=unique_orthogroups.index, y=unique_orthogroups.values, ax=axes)
    axes.set_title('Number of Unique and Shared Orthogroups', fontsize=16)
    axes.set_xlabel('Number of Species Sharing Orthogroups', fontsize=14)
    axes.set_ylabel('Count', fontsize=14)

    # Obtener el total de bins en el eje X
    total_bins = len(unique_orthogroups)

    # Establecer nbins como la mitad de los bins totales
    nbins = max(1, total_bins // 2)  # Aseguramos que al menos haya 1 bin

    # Limitar el número de etiquetas en el eje X a la mitad de los bins totales
    axes.xaxis.set_major_locator(MaxNLocator(integer=True, prune='both', nbins=nbins))

    # Guardar la imagen
    crear_directorio_plots()
    image_path = os.path.join('static', 'plots', 'figura_2.png')
    plt.savefig(image_path)
    plt.close(fig)

    return 'figura_2.png'

def generar_figura_3(gene_count_df, species, umbral=7):
    fig = plt.figure(figsize=(30, 17))  # Ajustamos el tamaño de la figura
    gs = GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[3, 1])  # Tamaño relativo de los subgráficos
    
    # Ejes
    ax1 = fig.add_subplot(gs[0, 0])  # Eje para el gráfico principal
    ax2 = fig.add_subplot(gs[0, 1])  # Eje para la leyenda de abreviaturas en formato tabla
    ax3 = fig.add_subplot(gs[1, :])  # Eje para la leyenda de combinaciones en formato tabla
    ax2.axis('off')
    ax3.axis('off')

    abreviaturas = {}
    abreviatura_usadas = set()

    for sp in species:
        abbr = ''.join([p[0].upper() for p in sp.split()[:2]])
        
        # Si la abreviatura ya ha sido usada, le agregamos un contador para hacerla única
        contador = 1
        abbr_unica = abbr
        while abbr_unica in abreviatura_usadas:
            abbr_unica = f"{abbr}{contador}"
            contador += 1
        
        abreviaturas[sp] = abbr_unica
        abreviatura_usadas.add(abbr_unica)

    gene_count_df = gene_count_df[(gene_count_df[species] > 0).any(axis=1)]
    gene_counts = {sp: 0 for sp in species}
    combinaciones_generadas = {}

    for _, row in gene_count_df.iterrows():
        presentes = row[species][row[species] > 0].index.tolist()
        if len(presentes) == 1:
            gene_counts[presentes[0]] += row[presentes[0]]
        elif len(presentes) > 1:
            combinacion = ' + '.join(sorted(presentes))
            if combinacion in combinaciones_generadas:
                combinaciones_generadas[combinacion] += row[presentes].sum()
            else:
                combinaciones_generadas[combinacion] = row[presentes].sum()

    gene_counts.update(combinaciones_generadas)
    total_genes = sum(gene_counts.values())
    porcentajes = {k: (v / total_genes) * 100 for k, v in gene_counts.items() if (v / total_genes) * 100 > 0.5}
    etiquetas_abreviadas = {}
    combinaciones_numeradas = {}
    contador_combinacion = 1

    for combinacion, porcentaje in porcentajes.items():
        especies = combinacion.split(' + ')
        if len(especies) > umbral:
            etiqueta_numerada = f'Combination {contador_combinacion}'
            combinaciones_numeradas[etiqueta_numerada] = ' + '.join([abreviaturas[sp] for sp in especies])
            etiquetas_abreviadas[combinacion] = etiqueta_numerada
            contador_combinacion += 1
        else:
            etiquetas_abreviadas[combinacion] = ' + '.join([abreviaturas[sp] for sp in especies])

    n_especies = {k: len(k.split(' + ')) for k in porcentajes.keys()}
    utilizados = sorted(set(n_especies.values()))

    ### Asignación de colores usando Viridis para combinaciones de especies
    cmap_combinaciones = colormaps.get_cmap('viridis')  # Sin el tercer argumento

    # Generar los colores utilizando un rango basado en el número de "utilizados"
    colores_combinaciones = [cmap_combinaciones(i / (len(utilizados) - 1)) for i in range(len(utilizados))]

    # Asignar los colores a cada combinación de especies
    colores = {}
    for i, n in enumerate(utilizados):
        color = colores_combinaciones[i]  # Asignar color de Viridis
        for k in n_especies:
            if n_especies[k] == n:
                colores[k] = color  # Asignar el color de Viridis a cada combinación de especies

    ### Asignación de colores claros manuales para la tabla de abreviaturas
    colores_claros = ["#add8e6", "#90ee90", "#ffb6c1", "#ffa07a", "#f0e68c", "#e0ffff", "#fafad2", "#d3d3d3", "#ffefd5", "#ffdab9",
                    "#e6e6fa", "#dda0dd", "#b0e0e6", "#bc8f8f", "#f5f5dc", "#ffe4e1", "#d8bfd8", "#d2b48c", "#add8e6", "#deb887"]
    color_cycle_claros = cycle(colores_claros)

    abreviaturas_cortas = {sp: sp[:5] for sp in species}  # Abreviaturas cortas (5 caracteres)
    unique_groups = list(set(abreviaturas_cortas.values()))  # Grupos únicos basados en las abreviaturas
    group_colors = {group: next(color_cycle_claros) for group in unique_groups}  # Asignamos un color claro a cada grupo

    ### Asignación de colores usando Viridis para los grupos de abreviaturas
    cmap_grupos = colormaps.get_cmap('viridis')  # Sin el tercer argumento
    colores_grupos = [cmap_grupos(i / (len(unique_groups) - 1)) for i in range(len(unique_groups))]


    # Crear la tabla de abreviaturas
    tabla_abreviaturas = pd.DataFrame(list(abreviaturas.items()), columns=['Strain Name', 'Abbreviation'])
    tabla = ax2.table(cellText=tabla_abreviaturas.values, cellLoc='center', loc='center')

    # Obtener la posición y tamaño de la tabla para ajustar la posición del título dinámicamente
    renderer = fig.canvas.get_renderer()
    tabla_pos = tabla.get_window_extent(renderer)

    # Calcular la nueva posición del título basado en la altura de la tabla
    # Queremos que esté justo por encima de la tabla
    y_title_position = tabla_pos.ymax / fig.bbox.ymax + 0.06 #Aquí ajusto bien para que el título Abbreviations se encuentra justo al lado de la tabla.

    # Colocar el título con anotación directamente sobre la tabla (centrado y ajustado en Y)
    ax2.annotate('Abbreviations', xy=(0.8, y_title_position), xycoords='figure fraction', 
                ha='center', fontsize=16, fontweight='bold', annotation_clip=False)

    # Aplicar los colores y ajustar el tamaño de las celdas de la tabla
    for (i, j), cell in tabla.get_celld().items():
        strain_name = tabla_abreviaturas.iloc[i, 0] if i < len(tabla_abreviaturas) else ""
        group_key = strain_name[:5]
        if j == 0:
            cell.set_fontsize(9)
        elif j == 1:
            cell.set_fontsize(14)
        if group_key in group_colors:
            cell.set_facecolor(group_colors[group_key])

    # Escalar la tabla para que el tamaño de la fuente y el espaciado sean adecuados
    tabla.scale(1.2, 1.2)

    # Eliminar esta línea que genera el segundo título redundante
    # ax2.set_title('Abbreviations', fontsize=14, pad=10)  # Eliminar esta línea
    ax2.axis('off')

    # Gráfico principal en ax1
    bars = ax1.barh([etiquetas_abreviadas[k] for k in porcentajes.keys()], list(porcentajes.values()), color=[colores[k] for k in porcentajes.keys()])
    ax1.set_xlabel('Percentage of Proteins (%)', fontweight='bold',)
    ax1.set_title('Percentage of Unique and Shared Proteins between Species', fontweight='bold',)

    for bar in bars:
        width = bar.get_width()
        ax1.text(width + 0.1, bar.get_y() + bar.get_height()/2, f'{width:.2f}%', va='center', rotation=10)

    # Personalizar los ejes: ocultar superior y derecho, pero mantener los ejes X e Y
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Eliminar duplicados en la leyenda y asignar correctamente los colores
    unique_legend_labels = []
    species_legend_patches = []
    legend_colors = {}
    for k, n in n_especies.items():
        if n not in unique_legend_labels:
            legend_colors[n] = colores[k]
            unique_legend_labels.append(n)

    # Ordenar la leyenda de menos a más especies
    sorted_legend_labels = sorted(unique_legend_labels)

    # Crear los parches de la leyenda ordenados
    species_legend_patches = [plt.Line2D([0], [0], marker='o', color='w', label=f'{n} species', markersize=10, markerfacecolor=legend_colors[n]) for n in sorted_legend_labels]
    ax1.legend(handles=species_legend_patches, title="Number of Species", loc='upper right')


    # Tabla de combinaciones (Figura 3) - Solo si hay combinaciones
    if combinaciones_numeradas:
        tabla_combinaciones = pd.DataFrame(list(combinaciones_numeradas.items()), columns=['Combination', 'Species'])
        tabla3 = ax3.table(cellText=tabla_combinaciones.values, cellLoc='center', loc='center')

        # Obtener la posición y tamaño de la tabla para ajustar la posición del título dinámicamente
        renderer = fig.canvas.get_renderer()
        tabla_pos = tabla3.get_window_extent(renderer)

        # Calcular la nueva posición del título basado en la altura de la tabla
        # Queremos que esté justo por encima de la tabla
        y_title_position = tabla_pos.ymax / fig.bbox.ymax - 0.005   # Ajustar el 0.02 para afinar la posición si es necesario

        # Colocar el título con anotación directamente sobre la tabla (centrado y ajustado en Y)
        ax3.annotate('Combinations', xy=(0.55, y_title_position), xycoords='figure fraction', 
                    ha='center', fontsize=16, fontweight='bold', annotation_clip=False)

        # Aplicar los colores que están en la Figura 1 (colores) a las filas correspondientes de la tabla
        for (i, j), cell in tabla3.get_celld().items():
            if i < len(tabla_combinaciones):  # Solo aplicar a las filas válidas
                # Obtener el número de especies en la combinación (segunda columna)
                species_combination = tabla_combinaciones.iloc[i, 1]
                num_species = len(species_combination.split(' + '))
                # Aplicar el color correspondiente al número de especies
                for k in colores:
                    if n_especies[k] == num_species:
                        cell.set_facecolor(colores[k])  # Asignar el color correspondiente
                # Ajustar el ancho de las columnas
                if j == 0:  # Primera columna (Combination)
                    cell.set_width(0.10)
                elif j == 1:  # Segunda columna (Species)
                    cell.set_width(0.85)

        # Escalar la tabla
        tabla3.scale(1.2, 1.2)

        ax3.axis('off')

    plt.tight_layout()

    crear_directorio_plots()
    
    image_path = os.path.join('static', 'plots', 'figura_3.png')
    plt.savefig(image_path)
    plt.close(fig)

    return 'figura_3.png', abreviaturas

def crear_upset_plot_ortogrupos(gene_count_df, species, title):
    """Genera un UpSet plot de ortogrupos basado en las especies seleccionadas y devuelve las combinaciones de ortogrupos."""
    memberships = []
    ortogrupos_por_combinacion = {}

    # Construir las combinaciones de ortogrupos presentes en cada especie
    for _, row in gene_count_df[species].iterrows():
        combinacion = tuple([sp for sp, val in zip(species, row) if val > 0])
        memberships.append(combinacion)

        # Contar cuántos ortogrupos pertenecen a cada combinación
        if combinacion in ortogrupos_por_combinacion:
            ortogrupos_por_combinacion[combinacion].append(row.name)
        else:
            ortogrupos_por_combinacion[combinacion] = [row.name]

    # Crear los datos para el UpSet plot
    upset_data = from_memberships(memberships)

    # Crear el gráfico UpSet
    fig = plt.figure(figsize=(12, 8))
    upset = UpSet(
        upset_data,
        subset_size='count',
        show_counts='%d',
        sort_by='degree',  # Ordenar las intersecciones por grado (número de especies compartidas)
    )

    # Dibujar el gráfico
    upset.plot(fig=fig)
    plt.suptitle(title, fontsize=16)

    # Guardar la imagen en el directorio de plots
    crear_directorio_plots()
    image_path = os.path.join('static', 'plots', f'{title.replace(" ", "_")}.png')
    plt.savefig(image_path, bbox_inches='tight')
    plt.close()

    return f'{title.replace(" ", "_")}.png', ortogrupos_por_combinacion

def crear_upset_plot_proteinas(gene_count_df, species, title, ortogrupos_por_combinacion):
    """Genera un UpSet plot de proteínas basado en las combinaciones de ortogrupos del primer plot."""
    memberships = []
    protein_counts = {}

    # Para cada combinación de ortogrupos, contar las proteínas presentes en las especies correspondientes
    for combinacion, ortogrupos in ortogrupos_por_combinacion.items():
        total_proteins = 0

        # Sumar las proteínas presentes en todos los ortogrupos de la combinación
        for ortogrupo in ortogrupos:
            protein_count_per_species = gene_count_df.loc[ortogrupo, species].sum()
            total_proteins += protein_count_per_species

        memberships.append(combinacion)
        protein_counts[combinacion] = total_proteins

    # Crear los datos para el UpSet plot basado en proteínas usando los conteos calculados
    upset_data = from_memberships(memberships, data=list(protein_counts.values()))

    # Crear el gráfico UpSet
    fig = plt.figure(figsize=(12, 8))
    upset = UpSet(
        upset_data,
        subset_size='sum',
        show_counts='%d',
        sort_by='degree',  # Ordenar las intersecciones por grado (número de especies compartidas)
    )

    # Dibujar el gráfico
    upset.plot(fig=fig)
    plt.suptitle(title, fontsize=16)

    # Guardar la imagen en el directorio de plots
    crear_directorio_plots()
    image_path = os.path.join('static', 'plots', f'{title.replace(" ", "_")}.png')
    plt.savefig(image_path, bbox_inches='tight')
    plt.close()

    return f'{title.replace(" ", "_")}.png', protein_counts

def generar_archivo_excel_upsetplots_v2(ortogrupos_por_combinacion, orthogroups_df, species, abreviaturas, filename='UpSetPlot_Data_v2.xlsx'):
    """Genera un archivo Excel con múltiples hojas, donde cada hoja representa una combinación de especies.
    Cada hoja contiene una lista de ortogrupos y las proteínas correspondientes para cada especie en columnas separadas."""

    with pd.ExcelWriter(filename) as writer:
        for combinacion, ortogrupos in ortogrupos_por_combinacion.items():
            # Crear una lista para almacenar los datos de ortogrupos y proteínas para cada especie
            combinacion_data = []

            # Iterar sobre los ortogrupos y extraer el Orthogroup ID y las proteínas para cada especie
            for ortogrupo in ortogrupos:
                row_data = {'Orthogroup ID': ortogrupo}

                # Agregar las proteínas para cada especie en las columnas correspondientes
                for sp in species:
                    if sp in orthogroups_df.columns:
                        proteins = orthogroups_df.at[ortogrupo, sp] if pd.notna(orthogroups_df.at[ortogrupo, sp]) else ""
                        row_data[sp] = proteins

                combinacion_data.append(row_data)

            # Crear un DataFrame para la hoja de la combinación actual
            df_combinacion = pd.DataFrame(combinacion_data)

            # Usar las abreviaturas para crear el nombre de la hoja
            hoja_nombre = ' + '.join([abreviaturas.get(sp, sp) for sp in combinacion])
            if not hoja_nombre:  # Verificar si el nombre de la hoja está vacío
                hoja_nombre = 'Unnamed_Combination'  # Asignar un nombre predeterminado en caso de que esté vacío
            hoja_nombre = hoja_nombre[:31]  # Limitar a 31 caracteres para cumplir con el límite de Excel

            # Escribir el DataFrame de la combinación en una nueva hoja del archivo Excel
            df_combinacion.to_excel(writer, sheet_name=hoja_nombre, index=False)

    return filename

def filter_orthogroups(orthogroups_df, uniprot_ids):
    """Filtra ortogrupos que contienen genes con IDs de UniProt proporcionados."""
    selected_orthogroups = []
    for idx, row in orthogroups_df.iterrows():
        genes = set()
        for col in orthogroups_df.columns[1:]:
            if pd.notna(row[col]):
                extracted_ids = [gene.split('|')[1] for gene in row[col].split(', ') if '|' in gene]
                genes.update(extracted_ids)
        if genes & set(uniprot_ids):
            selected_orthogroups.append(row['Orthogroup'])

    print(f'Ortogrupos seleccionados: {selected_orthogroups}')  # Añadir esto
    return selected_orthogroups

def graficar_upset_plots_proteome(orthogroups_df, selected_orthogroups, species):
    """Genera UpSet plots para ortogrupos y genes específicos del proteoma."""
    data = orthogroups_df[orthogroups_df['Orthogroup'].isin(selected_orthogroups)]

    gene_presence = pd.DataFrame(index=data['Orthogroup'])
    total_genes = data.set_index('Orthogroup')

    for sp in species:
        total_genes[sp] = total_genes[sp].str.split(', ').apply(lambda x: len(set(x)) if isinstance(x, list) else 0)
        gene_presence[sp] = total_genes[sp] > 0
    
    total_genes = total_genes.apply(pd.to_numeric, errors='coerce')
    gene_presence = gene_presence.astype(int)

    memberships = []
    gene_counts = total_genes.sum(axis=1)  # Sumar las proteínas para cada ortogrupo

    for index, row in gene_presence.iterrows():
        memberships.append([sp for sp in species if row[sp] > 0])

    # Crear el UpSet plot para ortogrupos
    upset_data_groups = from_memberships(memberships)
    fig_ortogrupos = plt.figure(figsize=(12, 8))
    upset = UpSet(upset_data_groups, subset_size='count', show_counts=True)
    upset.plot(fig=fig_ortogrupos)
    plt.suptitle('UpSet Plot of Orthogroups (Proteome)')
    image_path_ortogrupos = os.path.join('static', 'plots', 'upset_plot_ortogrupos.png')
    plt.savefig(image_path_ortogrupos, bbox_inches='tight')
    plt.close(fig_ortogrupos)

    # Crear el UpSet plot para genes
    upset_data_genes = from_memberships(memberships, data=gene_counts)
    fig_genes = plt.figure(figsize=(12, 8))
    upset = UpSet(upset_data_genes, subset_size='sum', show_counts=True)
    upset.plot(fig=fig_genes)
    plt.suptitle('UpSet Plot of Genes (Proteome)')
    image_path_genes = os.path.join('static', 'plots', 'upset_plot_genes.png')
    plt.savefig(image_path_genes, bbox_inches='tight')
    plt.close(fig_genes)

    # Devolver solo los nombres de los archivos, no las rutas completas
    return 'upset_plot_ortogrupos.png', 'upset_plot_genes.png'

def generar_archivo_excel_upsetplots_filtrado(orthogroups_df, ortogrupos_por_combinacion, species, abreviaturas, selected_orthogroups, filename='UpSetPlot_Filtrado.xlsx'):
    """Genera un archivo Excel filtrado con los selected_orthogroups proporcionados."""

    # Extraer la parte numérica de los selected_orthogroups
    selected_numeric_ids = [int(re.search(r'\d+', og).group()) for og in selected_orthogroups]

    # Generar el archivo Excel completo
    with pd.ExcelWriter('UpSetPlot_Data_v2.xlsx') as writer:
        for combinacion, ortogrupos in ortogrupos_por_combinacion.items():
            combinacion_data = []

            for ortogrupo in ortogrupos:
                row_data = {'Orthogroup ID': ortogrupo}

                for sp in species:
                    if sp in orthogroups_df.columns:
                        proteins = orthogroups_df.at[ortogrupo, sp] if pd.notna(orthogroups_df.at[ortogrupo, sp]) else ""
                        row_data[sp] = proteins

                combinacion_data.append(row_data)

            df_combinacion = pd.DataFrame(combinacion_data)

            # Usar las abreviaturas para crear el nombre de la hoja
            hoja_nombre = ' + '.join([abreviaturas.get(sp, sp) for sp in combinacion])

            # Si el nombre de la hoja está vacío, asignar un nombre predeterminado
            if not hoja_nombre:
                hoja_nombre = 'Unnamed_Combination'

            # Limitar a 31 caracteres para cumplir con el límite de Excel
            hoja_nombre = hoja_nombre[:31]

            # Escribir el DataFrame de la combinación en una nueva hoja del archivo Excel
            df_combinacion.to_excel(writer, sheet_name=hoja_nombre, index=False)

    # Cargar el archivo Excel completo y filtrar
    original_data = pd.read_excel('UpSetPlot_Data_v2.xlsx', sheet_name=None)

    # Diccionario para almacenar las hojas filtradas
    filtered_sheets = {}

    # Filtrar cada hoja del archivo original
    for sheet_name, df in original_data.items():
        # Filtrar las filas que contienen los ortogrupos seleccionados
        filtered_df = df[df['Orthogroup ID'].isin(selected_numeric_ids)]

        # Si la hoja no queda vacía, la añadimos al diccionario
        if not filtered_df.empty:
            filtered_sheets[sheet_name] = filtered_df

    # Crear el nuevo archivo Excel con solo las hojas filtradas
    with pd.ExcelWriter(filename) as writer:
        for sheet_name, df in filtered_sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    return filename

def read_orthogroups_data(zip_path):
    """Leer y devolver el contenido del archivo Orthogroups.tsv desde un archivo ZIP."""
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        with zip_ref.open('Orthogroups.tsv') as file:
            return pd.read_csv(file, sep='\t')

def generate_go_excel(zip_path, especies_interes, output_excel_path):
    # Utilizar la nueva función para leer datos
    ortho_df = read_orthogroups_data(zip_path)

    # Total de ortogrupos iniciales
    ortogrupos_iniciales = ortho_df.copy()

    # Calcular el porcentaje de anotación para cada ortogrupo en ortogrupos_iniciales
    porcentajes_anotacion = []
    for idx, row in ortogrupos_iniciales.iterrows():
        total_proteinas = row[1:].notna().sum()  # Número total de proteínas en el ortogrupo
        proteinas_interes = row[especies_interes].notna().sum()  # Número de proteínas de las especies de interés
        porcentaje_anotacion = (proteinas_interes / total_proteinas) * 100 if total_proteinas > 0 else 0
        porcentajes_anotacion.append(porcentaje_anotacion)

    # Añadir el porcentaje de anotación como una nueva columna en ortogrupos_iniciales
    ortogrupos_iniciales.insert(1, 'Annotation Percentage', porcentajes_anotacion)

    # Filtrar los ortogrupos que contienen al menos una proteína de alguna de las especies de interés
    ortogrupos_filtrados = ortho_df[ortho_df[especies_interes].notna().any(axis=1)].copy()
    ortogrupos_eliminados = ortho_df[~ortho_df.index.isin(ortogrupos_filtrados.index)].copy()
    ortogrupos_filtrados_interes = ortogrupos_filtrados[["Orthogroup"] + especies_interes].copy()
    ortogrupos_filtrados_interes['Annotation Percentage'] = [p for p in porcentajes_anotacion if p > 0]

    # Guardar los resultados en un archivo Excel
    with pd.ExcelWriter(output_excel_path) as writer:
        ortogrupos_iniciales.to_excel(writer, sheet_name='Initial Groups', index=False)
        ortogrupos_filtrados.to_excel(writer, sheet_name='Filtered Groups', index=False)
        ortogrupos_eliminados.to_excel(writer, sheet_name='Removed Groups', index=False)
        ortogrupos_filtrados_interes.to_excel(writer, sheet_name='Groups of Interest', index=False)

    return output_excel_path

def generate_go_image(zip_path, especies_interes, output_image_path):

    # Utilizar la nueva función para leer datos
    ortho_df = read_orthogroups_data(zip_path)

    # Calcular el porcentaje de anotación
    porcentajes_anotacion = []
    for idx, row in ortho_df.iterrows():
        total_proteinas = row[1:].notna().sum()
        proteinas_interes = row[especies_interes].notna().sum()
        porcentaje_anotacion = (proteinas_interes / total_proteinas) * 100 if total_proteinas > 0 else 0
        porcentajes_anotacion.append(porcentaje_anotacion)

    # Obtener los datos de anotación de ortogrupos
    all_annotation_percentages = pd.Series(porcentajes_anotacion)
    non_zero_annotation_percentages = all_annotation_percentages[all_annotation_percentages > 0]
    bins = range(0, 101, 10)  # De 0 a 100 con pasos de 10

    # Crear la cuadrícula de 2x2 para las gráficas
    fig, axes = plt.subplots(2, 2, figsize=(18, 16))
    sns.boxplot(y=all_annotation_percentages, color="skyblue", width=0.5, ax=axes[0, 0])
    median_all = all_annotation_percentages.median()
    q1_all = all_annotation_percentages.quantile(0.25)
    q3_all = all_annotation_percentages.quantile(0.75)
    axes[0, 0].axhline(median_all, color="orange", linestyle="--", label=f"Median: {median_all:.2f}%")
    axes[0, 0].axhline(q1_all, color="green", linestyle="--", label=f"25th Percentile (Q1): {q1_all:.2f}%")
    axes[0, 0].axhline(q3_all, color="purple", linestyle="--", label=f"75th Percentile (Q3): {q3_all:.2f}%")
    axes[0, 0].set_ylabel('Annotation Percentage')
    axes[0, 0].set_title('All Orthogroups')
    axes[0, 0].legend()
    axes[0, 0].text(-0.1, 1, 'a)', transform=axes[0, 0].transAxes, size=14, weight='bold')

    sns.boxplot(y=non_zero_annotation_percentages, color="lightgreen", width=0.5, ax=axes[0, 1])
    median_non_zero = non_zero_annotation_percentages.median()
    q1_non_zero = non_zero_annotation_percentages.quantile(0.25)
    q3_non_zero = non_zero_annotation_percentages.quantile(0.75)
    axes[0, 1].axhline(median_non_zero, color="orange", linestyle="--", label=f"Median: {median_non_zero:.2f}%")
    axes[0, 1].axhline(q1_non_zero, color="green", linestyle="--", label=f"25th Percentile (Q1): {q1_non_zero:.2f}%")
    axes[0, 1].axhline(q3_non_zero, color="purple", linestyle="--", label=f"75th Percentile (Q3): {q3_non_zero:.2f}%")
    axes[0, 1].set_ylabel('Annotation Percentage')
    axes[0, 1].set_title('Orthogroups with >0% Annotation')
    axes[0, 1].legend()
    axes[0, 1].text(-0.1, 1, 'b)', transform=axes[0, 1].transAxes, size=14, weight='bold')

    sns.histplot(all_annotation_percentages, bins=bins, color="cornflowerblue", kde=True, edgecolor="black", ax=axes[1, 0])
    axes[1, 0].set_xticks(bins)
    axes[1, 0].set_xlabel('Annotation Percentage')
    axes[1, 0].set_ylabel('Number of Orthogroups')
    axes[1, 0].set_title('Annotation Percentage Distribution - All Orthogroups')
    axes[1, 0].text(-0.1, 1, 'c)', transform=axes[1, 0].transAxes, size=14, weight='bold')

    for patch in axes[1, 0].patches:
        height = patch.get_height()
        if height > 0:
            axes[1, 0].text(patch.get_x() + patch.get_width() / 2, height + 50, f'{int(height)}', ha='center', va='bottom')

    sns.histplot(non_zero_annotation_percentages, bins=bins, color="lightgreen", kde=True, edgecolor="black", ax=axes[1, 1])
    axes[1, 1].set_xticks(bins)
    axes[1, 1].set_xlabel('Annotation Percentage')
    axes[1, 1].set_ylabel('Number of Orthogroups')
    axes[1, 1].set_title('Annotation Percentage Distribution - Orthogroups with >0% Annotation')
    axes[1, 1].text(-0.1, 1, 'd)', transform=axes[1, 1].transAxes, size=14, weight='bold')

    # Añadir etiquetas de frecuencia en el histograma de ortogrupos >0%
    for patch in axes[1, 1].patches:
        height = patch.get_height()
        if height > 0:
            axes[1, 1].text(
                patch.get_x() + patch.get_width() / 2,
                height + 50,
                f'{int(height)}',
                ha='center',
                va='bottom'
            )

    plt.tight_layout()
    fig.savefig(output_image_path)
    plt.close()

    return output_image_path

# Inicializar GODag después de asegurar que go_obo_file existe
go_obo_file = os.path.join('GOA_files', 'go-basic.obo')
if os.path.exists(go_obo_file):
    godag = GODag(go_obo_file)
else:
    raise FileNotFoundError("El archivo go-basic.obo no se encontró en la ruta especificada.")

def download_gaf_files(uniprot_ids):
    # Inicialización de id2gos y species_map
    id2gos = {uid: set() for uid in uniprot_ids}
    species_map = {}
    goa_files_dir = 'GOA_files'
    os.makedirs(goa_files_dir, exist_ok=True)  # Crear directorio si no existe

    # URLs de descarga de archivos GAF
    species_urls = {
        'Pseudomonas_aeruginosa_PAO1': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/36.P_aeruginosa_LMG_12228.goa',
        'Pseudomonas_putida_KT2440': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/109.P_putida_KT2440.goa',
        'Staphylococcus_aureus_PS47': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/22608.S_aureus_NCTC_8325.goa',
        'Mycobacteroides_abscessus_ATCC_19977': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/30878.M_abscessus.goa',
        'Mycobacterium_tuberculosis_H37Rv': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/30.M_tuberculosis_ATCC_25618.goa',
        'Mycobacterium_smegmatis_NCTC_8159': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/25827.M_smegmatis.goa',
        'Burkholderia_multivorans_ATCC_17616': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/31241.B_multivorans_Tohoku_University.goa',
        'Burkholderia_pseudomallei_K96243': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/20343.B_pseudomallei_K96243.goa',
        'Burkholderia_gladioli_BSR3': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/266524.B_gladioli_(strain_BSR3).goa',
        'Stenotrophomonas_maltophilia_K279A': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/30999.S_maltophilia_K279a.goa',
        'Haemophilus_influenzae_ATCC_51907': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/21.H_influenzae_ATCC_51907.goa',
        'Nocardia_farcinica_IFM_10152': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/20447.N_farcinica.goa',
        'Inquilinus_limosus': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/4109793.I_limosus.goa',
        'Acinetobacter_baumannii_ATCC_19606': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/283653.A_baumannii_ATCC_19606_=_CIP_7034.goa',
        'Klebsiella_pneumoniae_HS11286': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/269360.K_pneumoniae_subsp_pneumoniae_HS11286.goa',
        'Escherichia_coli_K12': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/18.E_coli_MG1655.goa',
        'Escherichia_coli_O157H7': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/20.E_coli_Sakai.goa',
        'Pandoraea_sputorum': 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/4128204.P_sputorum.goa',
    }

    # Descargar archivos GAF
    for species, url in species_urls.items():
        filename = os.path.join(goa_files_dir, Path(url).name)
        if not os.path.exists(filename):  # Descargar solo si no existe
            retries = 3
            for attempt in range(retries):
                try:
                    print(f"Downloading {species} from {url}... (Try {attempt + 1}/{retries})")
                    response = requests.get(url, stream=True, timeout=10)  # Timeout agregado
                    response.raise_for_status()  # Lanzar excepción si hay error en la solicitud
                    with open(filename, 'wb') as f:
                        shutil.copyfileobj(response.raw, f)
                    print(f"Downloaded {filename}")
                    break  # Salir del bucle si la descarga es exitosa
                except requests.exceptions.RequestException as e:
                    print(f"Error downloading {url}: {e}")
                    if attempt < retries - 1:
                        print("Retrying...")
                        time.sleep(5)  # Esperar antes de reintentar
                    else:
                        print("Failed to download after multiple attempts. Moving on.")

    # Descargar el archivo go-basic.obo si no existe
    go_obo_file = os.path.join(goa_files_dir, "go-basic.obo")
    if not os.path.exists(go_obo_file):
        try:
            print("Downloading GO OBO file...")
            response = requests.get("http://purl.obolibrary.org/obo/go/go-basic.obo", timeout=10)
            response.raise_for_status()
            with open(go_obo_file, 'wb') as f:
                f.write(response.content)
            print("Downloaded GO OBO file.")
        except requests.exceptions.RequestException as e:
            print(f"Error downloading GO OBO file: {e}")

    return id2gos, species_map, go_obo_file

def generate_foreground_analysis(uniprot_ids, use_orthogroups, ortogrupos_iniciales, ortogrupos_filtrados_interes):
    """Realiza el análisis de foreground basado en UniProt IDs y opción de ortogrupos."""
    
    # Validar que haya UniProt IDs proporcionados
    if not uniprot_ids:
        raise ValueError("No UniProt IDs provided for analysis.")

    # Descargar y obtener id2gos, species_map y go_obo_file
    id2gos, species_map, go_obo_file = download_gaf_files(uniprot_ids)

    # Convertir UniProt IDs a conjunto para operaciones rápidas
    uniprot_set = set(uniprot_ids)
    selected_orthogroups = set()
    protein_set = set()  # Usar un conjunto para evitar duplicados

    if use_orthogroups:
        for _, row in ortogrupos_iniciales.iterrows():
            orthogroup_id = row['Orthogroup']
            proteins = row[1:]
            for protein in proteins.dropna():
                if not isinstance(protein, str):
                    protein = str(protein)
                protein_ids = re.findall(r'\|([^|]+)\|', protein)
                if any(uid in uniprot_set for uid in protein_ids):
                    selected_orthogroups.add(orthogroup_id)
                    break
        selected_orthogroups_list = list(selected_orthogroups)
        for ortogroup_id in selected_orthogroups_list:
            filtered_rows = ortogrupos_filtrados_interes[ortogrupos_filtrados_interes['Orthogroup'] == ortogroup_id]
            for _, row in filtered_rows.iterrows():
                proteins = row.dropna().astype(str)
                for protein in proteins:
                    if protein != 'Porcentaje de Anotación':
                        protein_ids = re.findall(r'\|([^|]+)\|', protein)
                        protein_set.update(protein_ids)
    else:
        protein_set.update(uniprot_ids)
        for _, row in ortogrupos_iniciales.iterrows():
            orthogroup_id = row['Orthogroup']
            proteins = row[1:]
            for protein in proteins.dropna():
                if not isinstance(protein, str):
                    protein = str(protein)
                protein_ids = re.findall(r'\|([^|]+)\|', protein)
                if any(uid in uniprot_set for uid in protein_ids):
                    selected_orthogroups.add(orthogroup_id)
                    break

    final_protein_list = list(protein_set)
    return final_protein_list, selected_orthogroups_list

def load_background(file_path):
    with open(file_path, 'r') as f:
        ids = {line.strip().strip('"').strip(',') for line in f if line.strip()}
    return ids
##################################################################################################################
def ensure_gaf_file(gaf_file, species):
    if not os.path.exists(gaf_file):
        # Lógica de descarga (similar a la que ya tienes)
        try:
            # Descarga el archivo, descomprime si es necesario
            pass
        except requests.exceptions.RequestException as e:
            print(f"Error al descargar {gaf_file}: {e}")
    return gaf_file

def filter_by_depth(results, min_depth=2):
    return [r for r in results if godag.query_term(r.GO).depth >= min_depth]

def serialize_go_result(result):
    return {
        "GO": result.GO,
        "name": result.goterm.name,
        "p_fdr_bh": result.p_fdr_bh
    }


def get_species_gaf_files():
    return {
        'Pseudomonas_aeruginosa_PAO1': 'GOA_files/36.P_aeruginosa_LMG_12228.goa',
        'Pseudomonas_putida_KT2440': 'GOA_files/109.P_putida_KT2440.goa',
        'Staphylococcus_aureus_PS47': 'GOA_files/22608.S_aureus_NCTC_8325.goa',
        'Mycobacteroides_abscessus_ATCC_19977': 'GOA_files/30878.M_abscessus.goa',
        'Mycobacterium_tuberculosis_H37Rv': 'GOA_files/30.M_tuberculosis_ATCC_25618.goa',
        'Mycobacterium_smegmatis_NCTC_8159': 'GOA_files/25827.M_smegmatis.goa',
        'Burkholderia_multivorans_ATCC_17616': 'GOA_files/31241.B_multivorans_Tohoku_University.goa',
        'Burkholderia_pseudomallei_K96243': 'GOA_files/20343.B_pseudomallei_K96243.goa',
        'Burkholderia_gladioli_BSR3': 'GOA_files/266524.B_gladioli_(strain_BSR3).goa',
        'Stenotrophomonas_maltophilia_K279A': 'GOA_files/30999.S_maltophilia_K279a.goa',
        'Haemophilus_influenzae_ATCC_51907': 'GOA_files/21.H_influenzae_ATCC_51907.goa',
        'Nocardia_farcinica_IFM_10152': 'GOA_files/20447.N_farcinica.goa',
        'Inquilinus_limosus': 'GOA_files/4109793.I_limosus.goa',
        'Acinetobacter_baumannii_ATCC_19606': 'GOA_files/283653.A_baumannii_ATCC_19606_=_CIP_7034.goa',
        'Klebsiella_pneumoniae_HS11286': 'GOA_files/269360.K_pneumoniae_subsp_pneumoniae_HS11286.goa',
        'Escherichia_coli_K12': 'GOA_files/18.E_coli_MG1655.goa',
        'Escherichia_coli_O157H7': 'GOA_files/20.E_coli_Sakai.goa',
        'Pandoraea_sputorum': 'GOA_files/4128204.P_sputorum.goa',
    }

def generate_go_figure(bp_results, cc_results, mf_results):
    """Genera un gráfico con los resultados de análisis GO para BP, CC y MF."""

    # Validar que hay resultados para graficar
    if not bp_results and not cc_results and not mf_results:
        raise ValueError("No significant GO results provided for plotting.")

    # Preparar los datos para las gráficas con longitud máxima de nombres
    def prepare_plot_data(results, max_label_length=20):
        go_terms = [
            (r['name'][:max_label_length] + '...' if len(r['name']) > max_label_length else r['name'])
            for r in results
        ]
        p_values = [-np.log10(r['p_fdr_bh']) for r in results]
        return go_terms, p_values

    bp_terms, bp_pvalues = prepare_plot_data(bp_results) if bp_results else (["No terms"], [0])
    cc_terms, cc_pvalues = prepare_plot_data(cc_results) if cc_results else (["No terms"], [0])
    mf_terms, mf_pvalues = prepare_plot_data(mf_results) if mf_results else (["No terms"], [0])

    # Crear la figura y definir el diseño
    fig = plt.figure(figsize=(20, 15))
    gs = GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[0.25, 0.75])

    # Gráfico para Biological Process (BP) - Panel izquierdo
    ax_bp = fig.add_subplot(gs[:, 0])  # Ocupa ambas filas
    ax_bp.barh(bp_terms, bp_pvalues, color='lightblue')
    ax_bp.set_xlabel('-log10(p-value)')
    ax_bp.set_ylabel('GO Terms')
    ax_bp.set_title('GO Enrichment Analysis - Biological Process')

    # Gráfico para Cellular Component (CC) - Panel superior derecho
    ax_cc = fig.add_subplot(gs[0, 1])
    ax_cc.barh(cc_terms, cc_pvalues, color='lightgreen')
    ax_cc.set_xlabel('-log10(p-value)')
    ax_cc.set_ylabel('GO Terms')
    ax_cc.set_title('GO Enrichment Analysis - Cellular Component')

    # Gráfico para Molecular Function (MF) - Panel inferior derecho
    ax_mf = fig.add_subplot(gs[1, 1])
    ax_mf.barh(mf_terms, mf_pvalues, color='lightcoral')
    ax_mf.set_xlabel('-log10(p-value)')
    ax_mf.set_ylabel('GO Terms')
    ax_mf.set_title('GO Enrichment Analysis - Molecular Function')

    # Ajustar diseño para evitar superposición de elementos
    plt.tight_layout()

    # Comprobar si la carpeta de salida existe y guardar la figura
    output_dir = os.path.join('static', 'plots')
    os.makedirs(output_dir, exist_ok=True)
    figure_name = 'go_analysis_figure.png'
    output_path = os.path.join(output_dir, figure_name)
    plt.savefig(output_path)
    plt.close(fig)

    print(f"GO analysis figure saved at {output_path}.")  # Debugging print
    return figure_name  # Devolver solo el nombre del archivo
    pass

def create_go_excel_report(bp_results, cc_results, mf_results):
    """Genera un archivo Excel con 3 hojas: BP, CC y MF, con columnas de Procedimiento y -log10(p-value)."""
    def prepare_excel_data(results):
        """Prepara los datos para el Excel con el procedimiento y -log10(p-value)."""
        data = []
        for r in results:
            data.append({
                "Procedure": r.goterm.name,  # Cambiado de r["name"] a r.goterm.name
                "-log10(p-value)": -np.log10(r.p_fdr_bh) if r.p_fdr_bh > 0 else 0  # Cambiado de r["p_fdr_bh"] a r.p_fdr_bh
            })
        return pd.DataFrame(data)

    # Preparar los datos para cada hoja
    bp_df = prepare_excel_data(bp_results)
    cc_df = prepare_excel_data(cc_results)
    mf_df = prepare_excel_data(mf_results)

    # Crear el directorio si no existe
    output_dir = 'static/data_folder'
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'go_enrichment_report.xlsx')

    # Crear el archivo Excel
    with pd.ExcelWriter(output_path) as writer:
        bp_df.to_excel(writer, sheet_name='BP', index=False)
        cc_df.to_excel(writer, sheet_name='CC', index=False)
        mf_df.to_excel(writer, sheet_name='MF', index=False)

    print(f"Excel file generated at {output_path}")
    return output_path
    pass

##################################################################################################################
def open_browser():
    webbrowser.open("http://127.0.0.1:5000")
############################################################################################################
############################################################################################################
###############################   AQUI EMPIEZA LA PARTE DE @APP ROUTE     ##################################
############################################################################################################
############################################################################################################

# Rutas para la web
@app.route('/')
def index():
    return render_template('index.html')

@app.route('/reanalyze')
def reanalyze():
    """Redirige al punto en que se generaron las figuras 1 y 2"""
    # Si las figuras ya han sido generadas, redirigimos al usuario a esa vista
    if 'figuras_generadas' in session:
        plot_url_1 = session['plot_url_1']
        plot_url_2 = session['plot_url_2']
        ruta_carpeta = session.get('ruta_carpeta', '')
        species = list(session.get('species', []))  # Recuperar las especies del dataframe
        return render_template('resultado.html', 
                               plot_url_1=plot_url_1, 
                               plot_url_2=plot_url_2, 
                               species=species, 
                               ruta_carpeta=ruta_carpeta)
    # Si no hay datos, redirigimos al inicio
    return redirect(url_for('index'))

@app.route('/cargar_carpeta', methods=['POST', 'GET'])
def cargar_carpeta():
    global orthogroups_df, gene_count_df, species_selected, ortogrupos_por_combinacion, abreviaturas

    # Ruta del archivo ZIP preseleccionado
    ruta_archivo_preseleccionado = os.path.join('static', 'Orthogroups.zip')

    if request.method == 'POST':
        # Procesar el archivo subido por el usuario
        if 'file' in request.files:
            folder = request.files['file']
            if folder.filename == '':
                flash('No se seleccionó ningún archivo ZIP.')
                return redirect(url_for('index'))
            if not folder.filename.endswith('.zip'):
                flash('El formato de archivo es incorrecto. Sube un archivo ZIP.')
                return redirect(url_for('index'))

            # Procesar el archivo ZIP subido
            ruta_carpeta = os.path.join('static', 'data_folder')
            os.makedirs(ruta_carpeta, exist_ok=True)
            folder_path = os.path.join(ruta_carpeta, folder.filename)
            folder.save(folder_path)

            try:
                import zipfile
                with zipfile.ZipFile(folder_path, 'r') as zip_ref:
                    zip_ref.extractall(ruta_carpeta)

                gene_count_df, orthogroups_df, single_copy_df, unassigned_df = cargar_datos_carpeta(ruta_carpeta)

                # Generar las figuras
                plot_url_1 = generar_figura_1(gene_count_df)
                plot_url_2 = generar_figura_2(gene_count_df)

                # Guardar en la sesión que las figuras ya fueron generadas
                session['figuras_generadas'] = True
                session['plot_url_1'] = plot_url_1
                session['plot_url_2'] = plot_url_2
                session['ruta_carpeta'] = ruta_carpeta
                session['species'] = list(gene_count_df.columns[1:])  # Guardar las especies para reutilizarlas

                return render_template('resultado.html', 
                                       plot_url_1=plot_url_1, 
                                       plot_url_2=plot_url_2, 
                                       species=list(gene_count_df.columns[1:]), 
                                       ruta_carpeta=ruta_carpeta)
            except Exception as e:
                flash(f"Error al procesar el archivo subido: {str(e)}")
                print(f"Error: {str(e)}")  # Mensaje en consola
                return redirect(url_for('index'))

        elif 'especies' in request.form:
            species_selected = request.form.getlist('especies')
            ruta_carpeta = request.form['ruta_carpeta']

            if not species_selected:
                flash('Por favor selecciona al menos una especie.')
                return redirect(url_for('index'))

            # Cargar datos para las especies seleccionadas
            gene_count_df, orthogroups_df, _, _ = cargar_datos_carpeta(ruta_carpeta)

            # Generar Figura 3
            img_path_3, abreviaturas = generar_figura_3(gene_count_df, species_selected)

            # Generar UpSet plots (Figuras 4 y 5)
            img_path_4, ortogrupos_por_combinacion = crear_upset_plot_ortogrupos(gene_count_df, species_selected, 'UpSet Plot of All Orthogroups')
            img_path_5, protein_counts = crear_upset_plot_proteinas(gene_count_df, species_selected, 'UpSet Plot of All Proteins', ortogrupos_por_combinacion)

            return render_template('figuras_adicionales.html', 
                                   img_path_3=img_path_3,
                                   img_path_4=img_path_4,
                                   img_path_5=img_path_5)

    # Si se utiliza el archivo preseleccionado
    if request.method == 'GET' and 'filename' in request.args:
        print("Solicitud GET recibida con filename:", request.args.get('filename'))
        if request.args.get('filename') == 'Orthogroups.zip':
            print("Procesando archivo preseleccionado Orthogroups.zip")
            ruta_carpeta = os.path.join('static', 'data_folder')
            os.makedirs(ruta_carpeta, exist_ok=True)

            # Usar la misma lógica que el archivo subido por el usuario
            try:
                import zipfile
                with zipfile.ZipFile(ruta_archivo_preseleccionado, 'r') as zip_ref:
                    zip_ref.extractall(ruta_carpeta)

                gene_count_df, orthogroups_df, single_copy_df, unassigned_df = cargar_datos_carpeta(ruta_carpeta)

                # Generar las figuras
                plot_url_1 = generar_figura_1(gene_count_df)
                plot_url_2 = generar_figura_2(gene_count_df)

                # Guardar en la sesión que las figuras ya fueron generadas
                session['figuras_generadas'] = True
                session['plot_url_1'] = plot_url_1
                session['plot_url_2'] = plot_url_2
                session['ruta_carpeta'] = ruta_carpeta
                session['species'] = list(gene_count_df.columns[1:])

                return render_template('resultado.html', 
                                       plot_url_1=plot_url_1, 
                                       plot_url_2=plot_url_2, 
                                       species=list(gene_count_df.columns[1:]), 
                                       ruta_carpeta=ruta_carpeta)
            except Exception as e:
                flash(f"Error al procesar el archivo preseleccionado: {str(e)}")
                print(f"Error: {str(e)}")  # Mensaje en consola
                return redirect(url_for('index'))

    # Verificar si las figuras ya fueron generadas anteriormente
    if 'figuras_generadas' in session:
        plot_url_1 = session['plot_url_1']
        plot_url_2 = session['plot_url_2']
        return render_template('resultado.html', 
                               plot_url_1=plot_url_1, 
                               plot_url_2=plot_url_2, 
                               species=list(gene_count_df.columns[1:]),
                               ruta_carpeta=session.get('ruta_carpeta', ''))

    return redirect(url_for('index'))

@app.route('/create_excel')
def create_excel():
    """Generar el archivo Excel cuando se solicite mediante el botón."""
    global orthogroups_df, species_selected, gene_count_df, ortogrupos_por_combinacion, abreviaturas

    # Verifica si los datos están cargados; si no, los carga
    if orthogroups_df is None:
        gene_count_df, orthogroups_df, _, _ = cargar_datos_carpeta('static/data_folder')

    # Crear el archivo Excel
    excel_filename = generar_archivo_excel_upsetplots_v2(ortogrupos_por_combinacion, orthogroups_df, species_selected, abreviaturas)

    return 'Excel generado'

@app.route('/download_excel')
def download_excel():
    excel_path = 'UpSetPlot_Data_v2.xlsx'
    return send_file(excel_path, as_attachment=True)

@app.route('/download/<filename>')
def download_image(filename):
    return send_file(os.path.join('static', 'plots', filename), as_attachment=True)

@app.route('/generate_new_figures', methods=['POST'])
def generate_new_figures():
    global orthogroups_df, species_selected, ortogrupos_por_combinacion, selected_orthogroups  # Añadir selected_orthogroups como global

    # Obtener los UniProt IDs enviados desde el frontend
    data = request.get_json()
    uniprot_ids = data.get('data', '').split()

    # Filtrar los ortogrupos utilizando los UniProt IDs
    selected_orthogroups = filter_orthogroups(orthogroups_df, uniprot_ids)

    if selected_orthogroups:
        # Generar las figuras filtradas (Figuras 6 y 7)
        img_path_ortogrupos, img_path_genes = graficar_upset_plots_proteome(
            orthogroups_df, selected_orthogroups, species_selected
        )

        # Devolver las rutas de las figuras al frontend
        return jsonify({'img_path_6': img_path_ortogrupos, 'img_path_7': img_path_genes})
    else:
        return jsonify({'error': 'No se encontraron ortogrupos coincidentes con los IDs de UniProt proporcionados.'}), 400

@app.route('/create_excel_proteome_filtrado', methods=['POST'])
def create_excel_proteome_filtrado():
    global orthogroups_df, ortogrupos_por_combinacion, species_selected, abreviaturas, selected_orthogroups

    # Verifica que selected_orthogroups ya se ha generado
    if not selected_orthogroups:
        return jsonify({'error': 'No se han generado ortogrupos filtrados aún.'}), 400

    # Generar el archivo Excel utilizando los selected_orthogroups
    excel_filename = generar_archivo_excel_upsetplots_filtrado(orthogroups_df, ortogrupos_por_combinacion, species_selected, abreviaturas, selected_orthogroups)

    return send_file(excel_filename, as_attachment=True)

@app.route('/generate_go_excel', methods=['POST'])
def generate_go_excel_route():
    try:
        # Definir las especies de interés
        especies_interes = [
            'Pseudomonas aeruginosa PAO1',
            'Pseudomonas putida KT2440',
            'Staphylococcus aureus PS47',
            'Mycobacteroides abscessus ATCC 19977',
            'Mycobacterium tuberculosis H37Rv',
            'Mycobacterium smegmatis NCTC 8159',
            'Burkholderia multivorans ATCC_17616',
            'Burkholderia pseudomallei K96243',
            'Burkholderia gladioli BSR3',
            'Stenotrophomonas maltophilia K279A',
            'Haemophilus influenzae ATCC 51907',
            'Nocardia farcinica IFM 10152',
            'Inquilinus limosus',
            'Acinetobacter baumannii ATCC 19606',
            'Klebsiella pneumoniae HS11286',
            'E.Coli K12',
            'E.Coli O157H7',
            'Pandoraea sputorum NCTC13161'
        ]
        zip_path = os.path.join('static', 'data_folder', 'Orthogroups.zip')
        output_excel_path = os.path.join('static', 'data_folder', 'Gene_Ontology_Analysis.xlsx')
        
        # Llamada a la función para generar el Excel
        excel_file_path = generate_go_excel(zip_path, especies_interes, output_excel_path)
        return send_file(excel_file_path, as_attachment=True)

    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/generate_go_image', methods=['POST'])
def generate_go_image_route():
    try:
        # Definir las especies de interés
        especies_interes = [
            'Pseudomonas aeruginosa PAO1',
            'Pseudomonas putida KT2440',
            'Staphylococcus aureus PS47',
            'Mycobacteroides abscessus ATCC 19977',
            'Mycobacterium tuberculosis H37Rv',
            'Mycobacterium smegmatis NCTC 8159',
            'Burkholderia multivorans ATCC_17616',
            'Burkholderia pseudomallei K96243',
            'Burkholderia gladioli BSR3',
            'Stenotrophomonas maltophilia K279A',
            'Haemophilus influenzae ATCC 51907',
            'Nocardia farcinica IFM 10152',
            'Inquilinus limosus',
            'Acinetobacter baumannii ATCC 19606',
            'Klebsiella pneumoniae HS11286',
            'E.Coli K12',
            'E.Coli O157H7',
            'Pandoraea sputorum NCTC13161'
        ]
        zip_path = os.path.join('static', 'data_folder', 'Orthogroups.zip')
        output_image_path = os.path.join('static', 'plots', 'Gene_Ontology_Annotation.png')
        
        # Llamada a la función para generar la imagen
        image_file_path = generate_go_image(zip_path, especies_interes, output_image_path)
        return jsonify({"image_file_path": os.path.basename(image_file_path)})

    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/foreground_analysis', methods=['POST'])
def foreground_analysis():
    """Endpoint para manejar el análisis de foreground basado en UniProt IDs."""
    try:
        # Obtener datos enviados por el cliente
        uniprot_ids = request.json.get('uniprot_ids', [])
        use_orthogroups = request.json.get('use_orthogroups', False)  # Valor del checkbox
        print("Received UniProt IDs:", uniprot_ids)  # Debugging print
        print("Use Orthogroups flag:", use_orthogroups)  # Debugging print

        # Validar entrada
        if not uniprot_ids:
            return jsonify({"error": "No UniProt IDs provided"}), 400

        # Verificar si el archivo existe antes de cargarlo
        excel_path = 'static/data_folder/Gene_Ontology_Analysis.xlsx'
        if not os.path.exists(excel_path):
            print("Gene Ontology Analysis file not found at:", excel_path)  # Debugging print
            return jsonify({"error": "Gene Ontology Analysis file not found."}), 404

        # Cargar las hojas necesarias del archivo
        try:
            print("Loading sheets from Excel file...")  # Debugging print
            ortogrupos_iniciales = pd.read_excel(excel_path, sheet_name='Initial Groups')
            ortogrupos_filtrados_interes = pd.read_excel(excel_path, sheet_name='Groups of Interest')
            print("Sheets loaded successfully.")  # Debugging print
        except Exception as e:
            print("Error reading Excel sheets:", str(e))  # Debugging print
            return jsonify({"error": f"Error reading Excel sheets: {str(e)}"}), 500

        # Ejecutar análisis de foreground
        print("Running foreground analysis...")  # Debugging print
        foreground_proteins, selected_orthogroups_list = generate_foreground_analysis(
            uniprot_ids=uniprot_ids,
            use_orthogroups=use_orthogroups,
            ortogrupos_iniciales=ortogrupos_iniciales,
            ortogrupos_filtrados_interes=ortogrupos_filtrados_interes
        )
        print("Foreground analysis completed.")  # Debugging print

        # Guardar los resultados en la sesión
        session['foreground_proteins'] = foreground_proteins
        session['selected_orthogroups'] = selected_orthogroups_list
        print("Foreground proteins and selected orthogroups saved to session.")  # Debugging print

        # Mostrar resultados en la terminal
        print("Selected Orthogroups List:", selected_orthogroups_list, flush=True)
        print("Final Protein List (uniprot_ids):", foreground_proteins, flush=True)

        # Devolver resultado en formato JSON
        return jsonify({
            "foreground_proteins": foreground_proteins,
            "message": "Foreground analysis completed successfully"
        })

    except ValueError as ve:
        print("ValueError encountered:", str(ve))  # Debugging print
        return jsonify({"error": str(ve)}), 400
    except Exception as e:
        print("Unexpected error:", str(e))  # Debugging print
        return jsonify({"error": f"An unexpected error occurred: {str(e)}"}), 500

@app.route('/background_analysis', methods=['POST'])
def background_analysis():
    """Realiza el análisis de background basado en la selección de opciones y manejo de UniProt IDs."""
    try:
        # Recibir datos del cliente (opciones de selección y UniProt IDs personalizados)
        data = request.json
        background_choice = data.get('background_choice')  # Opción seleccionada (1, 2, 3, 4 o 5)
        use_orthogroups = data.get('use_orthogroups', False)  # Booleano para ortogrupos
        custom_uniprot_ids = data.get('custom_uniprot_ids', [])

        # Cargar datos de ortogrupos iniciales y filtrados
        excel_path = os.path.join('static', 'data_folder', 'Gene_Ontology_Analysis.xlsx')
        ortogrupos_iniciales = pd.read_excel(excel_path, sheet_name='Initial Groups')
        ortogrupos_filtrados_interes = pd.read_excel(excel_path, sheet_name='Groups of Interest')

        # Definir la ruta base para los archivos de background
        background_dir = os.path.join('static', 'data_folder')
        background_files = {
            'mycobacterium': os.path.join(background_dir, "Background_Mycobacterium_formatted.txt"),
            'burkholderia': os.path.join(background_dir, "Background_Burkholderia_formatted.txt"),
            'pseudomonas': os.path.join(background_dir, "Background_Pseudomonas_formatted.txt")
        }

        # Directorio donde están almacenados los archivos GAF ya descargados
        gaf_dir = os.path.join('static', 'GOA_files')

        # Inicializar el conjunto de IDs de background
        background_ids = set()

        # Procesar la opción seleccionada para cargar los UniProt IDs del background
        if background_choice == '1':
            background_ids = load_background(background_files['mycobacterium'])
        elif background_choice == '2':
            background_ids = load_background(background_files['burkholderia'])
        elif background_choice == '3':
            background_ids = load_background(background_files['pseudomonas'])
        elif background_choice == '4' and custom_uniprot_ids:
            background_ids = set(custom_uniprot_ids)
        elif background_choice == '5':
            # Usar archivos GAF previamente descargados
            for file_name in os.listdir(gaf_dir):
                if file_name.endswith('.gaf.gz'):
                    gaf_file_path = os.path.join(gaf_dir, file_name)
                    with gzip.open(gaf_file_path, 'rt') as f:
                        gaf_reader = GafReader(f)
                        for protein_id, go_ids in gaf_reader.get_id2gos_nss().items():
                            background_ids.add(protein_id)
        else:
            return jsonify({"error": "Invalid background choice or missing data"}), 400

        # Opcional: Manejar ortogrupos si está habilitado
        if use_orthogroups and background_choice in ['1', '2', '3']:
            selected_orthogroups_bg = set()
            for _, row in ortogrupos_iniciales.iterrows():
                orthogroup_id = row['Orthogroup']
                proteins = row[1:]  # Todas las columnas de proteínas después de 'Orthogroup'
                for protein in proteins.dropna():
                    if not isinstance(protein, str):
                        protein = str(protein)
                    protein_ids = re.findall(r'\|([^|]+)\|', protein)
                    if any(bg_id in background_ids for bg_id in protein_ids):
                        selected_orthogroups_bg.add(orthogroup_id)
                        break

            # Expandir los UniProt IDs del background usando los ortogrupos seleccionados
            background_expanded_ids = set()
            for ortogroup_id in selected_orthogroups_bg:
                filtered_rows = ortogrupos_filtrados_interes[
                    ortogrupos_filtrados_interes['Orthogroup'] == ortogroup_id]
                for _, row in filtered_rows.iterrows():
                    proteins = row.dropna().astype(str)
                    for protein in proteins:
                        if protein != 'Porcentaje de Anotación':
                            protein_ids = re.findall(r'\|([^|]+)\|', protein)
                            background_expanded_ids.update(protein_ids)

            # Actualizar el conjunto de IDs del background con los UniProt IDs ampliados
            background_ids = background_expanded_ids

        # Mostrar el número de proteínas seleccionadas como background
        print(f"Background selected with {len(background_ids)} proteins after orthogroup expansion (if applicable).", flush=True)

        # Guardar el resultado de background en la sesión
        session['background_ids'] = list(background_ids)
        print("Background IDs saved to session.")  # Debugging print

        # Retornar los resultados como JSON
        return jsonify({
            "background_ids_count": len(background_ids),
            "background_ids": list(background_ids)
        })

    except Exception as e:
        print(f"Error during background analysis: {e}")  # Debugging print
        return jsonify({"error": str(e)}), 500

@app.route('/gene_ontology_analysis', methods=['POST'])
def gene_ontology_analysis():
    try:
        # Obtener 'uniprot_ids' y 'background_ids' desde la solicitud
        uniprot_ids = request.json.get('uniprot_ids', [])
        background_ids = request.json.get('background_ids', [])

        # Obtener valores de Depth, P.value y Number of Terms desde la solicitud o establecer valores predeterminados
        min_depth = request.json.get('min_depth', 2)  # Valor predeterminado: 2
        p_value_threshold = request.json.get('p_value', 0.05)  # Valor predeterminado: 0.05
        max_terms = request.json.get('max_terms', None)  # Valor predeterminado: None (mostrar todos)

        if not uniprot_ids:
            print("No foreground (uniprot_ids) provided for Gene Ontology Analysis.")
            return jsonify({"error": "No foreground (uniprot_ids) provided"}), 400

        if not background_ids:
            print("No background (background_ids) provided for Gene Ontology Analysis.")
            return jsonify({"error": "No background (background_ids) provided"}), 400

        print("Foreground UniProt IDs:", uniprot_ids)
        print("Background IDs:", background_ids)
        print(f"Received parameters - min_depth: {min_depth}, p_value_threshold: {p_value_threshold}, max_terms: {max_terms}")

        # Inicializar id2gos
        id2gos = {uid: set() for uid in uniprot_ids}
        species_map = {}

        # Cargar archivos GAF y mapear anotaciones GO
        species_gaf_files = get_species_gaf_files()
        print("Loading annotations from GAF files...")
        for species, gaf_file in species_gaf_files.items():
            print(f"Processing GAF file for species: {species}")
            gaf_reader = GafReader(gaf_file)
            ns2assc = gaf_reader.get_ns2assc()

            for namespace, associations in ns2assc.items():
                for protein_id, go_ids in associations.items():
                    if protein_id in background_ids:
                        if protein_id not in id2gos:
                            id2gos[protein_id] = set(go_ids)
                            species_map[protein_id] = species
                        else:
                            id2gos[protein_id].update(go_ids)

        print("GO annotations loaded. Running enrichment analysis...")

        # Ejecutar análisis de enriquecimiento de GO
        goea = GOEnrichmentStudy(
            background_ids,
            id2gos,
            godag,  # Asegúrate de que 'godag' esté correctamente cargado
            methods=["fdr_bh"],
            log=None
        )
        goea_results = goea.run_study(uniprot_ids)

        # Filtrar resultados significativos según el p-value proporcionado
        significant_results = [r for r in goea_results if r.p_fdr_bh < p_value_threshold]
        print(f"Found {len(significant_results)} significant GO terms after p-value filtering.")

        # Filtrar resultados por tipo de categoría GO y aplicar filtro de profundidad
        bp_results = filter_by_depth(
            [r for r in significant_results if r.goterm.namespace == 'biological_process'], min_depth
        )
        cc_results = filter_by_depth(
            [r for r in significant_results if r.goterm.namespace == 'cellular_component'], min_depth
        )
        mf_results = filter_by_depth(
            [r for r in significant_results if r.goterm.namespace == 'molecular_function'], min_depth
        )

        # Limitar el número de términos si max_terms está definido
        if max_terms:
            bp_results = bp_results[:max_terms]
            cc_results = cc_results[:max_terms]
            mf_results = mf_results[:max_terms]

        print(f"BP results: {len(bp_results)}")
        print(f"CC results: {len(cc_results)}")
        print(f"MF results: {len(mf_results)}")

        # Generar la figura de análisis GO basada en los resultados obtenidos
        figure_path = generate_go_figure(
            [serialize_go_result(r) for r in bp_results],
            [serialize_go_result(r) for r in cc_results],
            [serialize_go_result(r) for r in mf_results]
        )

        print("Figure generated successfully at:", figure_path)

        # Generar el archivo Excel con los resultados proporcionados
        excel_path = create_go_excel_report(bp_results, cc_results, mf_results)
        print("Excel file generated successfully at:", excel_path)

        return jsonify({
            "bp_results": [serialize_go_result(r) for r in bp_results],
            "cc_results": [serialize_go_result(r) for r in cc_results],
            "mf_results": [serialize_go_result(r) for r in mf_results],
            "message": "Gene Ontology Analysis completed successfully",
            "image_file_path": figure_path,  # Devolver la ruta de la imagen
            "excel_file_path": excel_path    # Devolver la ruta del archivo Excel
        })

    except Exception as e:
        print(f"Error during Gene Ontology Analysis: {e}")
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    threading.Timer(1.25, open_browser).start()  # Se abrirá el navegador automáticamente
    app.run(debug=True, use_reloader=False)
