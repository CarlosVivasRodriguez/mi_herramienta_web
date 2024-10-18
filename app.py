import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from flask import Flask, render_template, request, redirect, url_for, flash, send_file, session, jsonify  # Mueve jsonify aquí
from matplotlib.gridspec import GridSpec  # Para las subgráficas en figura 3
from itertools import cycle  # Para generar combinaciones de colores
from upsetplot import UpSet, from_memberships  # Para generar UpSet plots
import warnings  # Para ignorar advertencias que no necesitas en el contexto de la app
import io  # Para manejo de buffers en memoria
import matplotlib.colors as mcolors
import webbrowser  # Importa el módulo para abrir el navegador automáticamente
import threading  # Para manejar el temporizador que abre el navegador con retras

# Importamos y configuramos matplotlib para evitar errores de GUI al generar figuras
import matplotlib
matplotlib.use('Agg')  # Establece el backend 'Agg' para evitar errores de GUI en servidores

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

# Función para cargar los datos
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

# Función para eliminar archivos temporales
def eliminar_archivo(ruta_archivo):
    if os.path.exists(ruta_archivo):
        os.remove(ruta_archivo)

# Función para asegurarse de que el directorio 'plots' exista
def crear_directorio_plots():
    if not os.path.exists(os.path.join('static', 'plots')):
        os.makedirs(os.path.join('static', 'plots'))

# Función para generar la primera figura (barras amarillas)
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

# Función para generar la segunda figura (barras azules)
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

    crear_directorio_plots()
    
    image_path = os.path.join('static', 'plots', 'figura_2.png')
    plt.savefig(image_path)
    plt.close(fig)

    return 'figura_2.png'

# Función para generar la tercera figura con tres subgráficas y combinaciones
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

    # Paleta de colores de Vega (category30 + un color extra)
    vega_palette = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#66c2a5", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
                    "#9a6324", "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080",
                    "#ffffff", "#000000", "#d2f53c", "#fabebe", "#008080", "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#aaffc3"]
    color_cycle = cycle(vega_palette)

    # Asignar colores basados en el número de especies que forman cada combinación
    colores = {}
    especies_por_n = sorted(utilizados)  # Aseguramos el orden para asignar colores
    for i, n in enumerate(especies_por_n):
        color = next(color_cycle)
        for k in n_especies:
            if n_especies[k] == n:
                colores[k] = color

    # Figura 2: Colorear en base a las primeras 5 letras de las especies (con colores claros)
    colores_claros = ["#add8e6", "#90ee90", "#ffb6c1", "#ffa07a", "#f0e68c", "#e0ffff", "#fafad2", "#d3d3d3", "#ffefd5", "#ffdab9",
                      "#e6e6fa", "#dda0dd", "#b0e0e6", "#bc8f8f", "#f5f5dc", "#ffe4e1", "#d8bfd8", "#d2b48c", "#add8e6", "#deb887"]
    color_cycle_claros = cycle(colores_claros)

    abreviaturas_cortas = {sp: sp[:5] for sp in species}
    unique_groups = list(set(abreviaturas_cortas.values()))
    group_colors = {group: next(color_cycle_claros) for group in unique_groups}  # Asignamos un color claro a cada grupo

   
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

# Función para generar el archivo Excel
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

import re  # Importamos la librería para manejar expresiones regulares

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




def open_browser():
    webbrowser.open("http://127.0.0.1:5000")
#####################################################################
#####################################################################
#            AQUI EMPIEZA LA PARTE DE @APP ROUTE
#####################################################################
#####################################################################

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






if __name__ == '__main__':
    # Solo abre el navegador en el proceso principal (evita que se abra en recargas automáticas de Flask)
    threading.Timer(1.25, open_browser).start()  # Abre el navegador después de un pequeño retraso
    app.run(debug=True, use_reloader=False)  # Desactiva el recargador automático