# streamlit_app.py

import streamlit as st
import geopandas as gpd
import pandas as pd
import numpy as np
import pydeck as pdk
import tempfile
import zipfile
import os
import io
from shapely.geometry import MultiPoint, Polygon
import h3
from sklearn.neighbors import NearestNeighbors
from shapely.geometry import Point
from sklearn.neighbors import NearestNeighbors
from shapely.geometry import Point
import zipfile
import io
import os
from shapely.geometry import box





# Función para calcular moda
def moda(x):
    return x.mode().iloc[0] if not x.mode().empty else np.nan

# Función para generar ZIP on-demand
def generar_zip(gdf_clean):
    with tempfile.TemporaryDirectory() as tmpzipdir:
        out_shp = os.path.join(tmpzipdir, "rendimiento_limpio.shp")
        gdf_clean.to_file(out_shp)
        
        zip_path = os.path.join(tmpzipdir, "rendimiento_limpio.zip")
        with zipfile.ZipFile(zip_path, "w") as zipf:
            for ext in [".shp", ".dbf", ".shx", ".prj"]:
                fpath = out_shp.replace(".shp", ext)
                if os.path.exists(fpath):
                    zipf.write(fpath, arcname=os.path.basename(fpath))
        
        with open(zip_path, "rb") as f:
            return io.BytesIO(f.read())

# Configuración inicial
st.set_page_config(page_title="Mapa de rendimiento + H3", layout="wide")
st.title("Mapa de rendimiento + Tratamiento + Ambiente (Grilla H3)")

# Tabs
tabs = st.tabs(["1️⃣ Limpieza mapa de rendimiento", "2️⃣ Fusión H3 + Tratamiento + Ambiente + Rinde","3️⃣ Limpieza de borde por cambio de tratamiento (buffer topológico)"])

# Estado
if "gdf_clean" not in st.session_state:
    st.session_state["gdf_clean"] = None
    st.session_state["col_rinde"] = None
if "gdf_h3_final" not in st.session_state:
    st.session_state["gdf_h3_final"] = None

# ------------------ TAB 1 ------------------

with tabs[0]:
    st.header("1️⃣ Limpieza mapa de rendimiento")

    uploaded_files = st.file_uploader(
        "Cargar archivos SHAPE (.shp, .dbf, .shx, .prj)", 
        type=["shp", "dbf", "shx", "prj"],
        accept_multiple_files=True
    )

    if uploaded_files:
        with tempfile.TemporaryDirectory() as tmpdir:
            for file in uploaded_files:
                with open(os.path.join(tmpdir, file.name), "wb") as f:
                    f.write(file.getbuffer())
            
            shp_files = [f.name for f in uploaded_files if f.name.endswith(".shp")]
            if not shp_files:
                st.error("No se subió un archivo .shp.")
            else:
                shp_path = os.path.join(tmpdir, shp_files[0])
                gdf = gpd.read_file(shp_path)
                
                if gdf.geom_type.iloc[0] == "MultiPoint":
                    points = []
                    for geom in gdf.geometry:
                        if isinstance(geom, MultiPoint):
                            for point in geom.geoms:
                                new_row = gdf.loc[[0]].copy()
                                new_row.geometry = point
                                points.append(new_row)
                    gdf = pd.concat(points).reset_index(drop=True)
                
                st.success(f"Archivo cargado con {len(gdf)} registros.")

                col_rinde = st.selectbox("Seleccionar variable de rendimiento", gdf.columns)
                col_ancho = st.selectbox("Seleccionar variable de ancho de trabajo", gdf.columns)

                if st.button("Aplicar limpieza"):
                    rinde_var = gdf[col_rinde]
                    ancho_var = gdf[col_ancho]
                    
                    IQR_val = rinde_var.quantile(0.75) - rinde_var.quantile(0.25)
                    Q1 = rinde_var.quantile(0.25)
                    Q3 = rinde_var.quantile(0.75)
                    
                    lower_bound_rinde = Q1 - 2 * IQR_val
                    upper_bound_rinde = Q3 + 2 * IQR_val
                    
                    moda_ancho = moda(ancho_var)
                    min_ancho = 0.8 * moda_ancho
                    
                    gdf_clean = gdf[
                        (rinde_var >= lower_bound_rinde) &
                        (rinde_var <= upper_bound_rinde) &
                        (ancho_var >= min_ancho)
                    ].copy()
                    
                    gdf_clean = gdf_clean.to_crs(epsg=4326)

                    st.session_state["gdf_clean"] = gdf_clean
                    st.session_state["col_rinde"] = col_rinde

                    if gdf_clean.empty:
                        st.warning("⚠️ El shapefile resultante está vacío.")
                    else:
                        st.info(f"Registros luego de limpieza: {len(gdf_clean)}")
                        st.write(f"Geometrías en GDF limpio: {gdf_clean.geometry.geom_type.value_counts().to_dict()}")

    # Mostrar mapa con PyDeck
    if st.session_state["gdf_clean"] is not None and not st.session_state["gdf_clean"].empty:
        gdf_clean = st.session_state["gdf_clean"]
        col_rinde = st.session_state["col_rinde"]

        gdf_valid = gdf_clean[gdf_clean.geometry.notnull()].copy()
        gdf_valid["lon"] = gdf_valid.geometry.x
        gdf_valid["lat"] = gdf_valid.geometry.y
        gdf_valid["rango"] = pd.qcut(gdf_valid[col_rinde], 5, labels=False)

        palette = [
            [165, 0, 38],
            [254, 224, 139],
            [255, 255, 191],
            [166, 217, 106],
            [26, 150, 65]
        ]

        gdf_valid["color"] = gdf_valid["rango"].apply(lambda x: palette[int(x)])

        scatter_layer = pdk.Layer(
            "ScatterplotLayer",
            data=gdf_valid,
            get_position='[lon, lat]',
            get_fill_color='color',
            get_radius=10,
            pickable=True
        )

        tile_layer = pdk.Layer(
            "TileLayer",
            data=None,
            min_zoom=0,
            max_zoom=20,
            tile_size=256,
            url="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
        )

        view_state = pdk.ViewState(
            latitude=gdf_valid["lat"].mean(),
            longitude=gdf_valid["lon"].mean(),
            zoom=12,
            pitch=0
        )

        deck = pdk.Deck(
            layers=[tile_layer, scatter_layer],
            initial_view_state=view_state,
            map_style='mapbox://styles/mapbox/satellite-v9'
        )

        st.pydeck_chart(deck)

        if st.checkbox("Mostrar dataframe limpio"):
            st.dataframe(gdf_clean.drop(columns="geometry"))

        zip_buffer = generar_zip(gdf_clean)
        st.download_button(
            label="Descargar shapefile limpio",
            data=zip_buffer,
            file_name="rendimiento_limpio.zip",
            mime="application/zip"
        )

from sklearn.neighbors import NearestNeighbors
from shapely.geometry import Point

from sklearn.neighbors import NearestNeighbors
from shapely.geometry import Point
import zipfile
import io
import os

from sklearn.neighbors import NearestNeighbors
from shapely.geometry import Point
import zipfile
import io
import os



# TAB 2 - Fusión Grilla (subida) + Tratamiento + Ambiente + Rinde
from sklearn.neighbors import NearestNeighbors
import zipfile
import io
import os

with tabs[1]:
    st.header("2️⃣ Fusión H3 + Tratamiento + Ambiente + Rinde")

    uploaded_grid = st.file_uploader("Archivos SHP (grilla H3)", type=["shp", "dbf", "shx", "prj"], accept_multiple_files=True, key="grilla_h3")
    uploaded_tratamientos = st.file_uploader("Archivos SHP (tratamientos)", type=["shp", "dbf", "shx", "prj"], accept_multiple_files=True, key="tratamiento")
    uploaded_ambientes = st.file_uploader("Archivos SHP (ambientes)", type=["shp", "dbf", "shx", "prj"], accept_multiple_files=True, key="ambientes")

    if uploaded_grid and uploaded_tratamientos and uploaded_ambientes:
        with tempfile.TemporaryDirectory() as tmpdir:

            # --- Grilla
            for file in uploaded_grid:
                with open(os.path.join(tmpdir, file.name), "wb") as f:
                    f.write(file.getbuffer())
            shp_files = [f.name for f in uploaded_grid if f.name.endswith(".shp")]
            shp_path = os.path.join(tmpdir, shp_files[0])
            gdf_grid = gpd.read_file(shp_path).to_crs(epsg=4326)

            # --- Tratamientos
            for file in uploaded_tratamientos:
                with open(os.path.join(tmpdir, file.name), "wb") as f:
                    f.write(file.getbuffer())
            shp_files = [f.name for f in uploaded_tratamientos if f.name.endswith(".shp")]
            shp_path = os.path.join(tmpdir, shp_files[0])
            gdf_trat = gpd.read_file(shp_path).to_crs(epsg=4326)

            # --- Ambientes
            for file in uploaded_ambientes:
                with open(os.path.join(tmpdir, file.name), "wb") as f:
                    f.write(file.getbuffer())
            shp_files = [f.name for f in uploaded_ambientes if f.name.endswith(".shp")]
            shp_path = os.path.join(tmpdir, shp_files[0])
            gdf_amb = gpd.read_file(shp_path).to_crs(epsg=4326)

        # --- Selección de columnas
        col_tratamiento = st.selectbox("Seleccionar columna de tratamiento", gdf_trat.columns)
        col_amb = st.selectbox("Seleccionar columna de ambiente", gdf_amb.columns)

        # --- Agrupar ambientes
        ambientes_unicos = sorted(gdf_amb[col_amb].dropna().unique())
        ambiente_grupos = {}

        st.write("Asignar cada ambiente a un grupo (1, 2 o 3):")
        for amb in ambientes_unicos:
            grupo = st.selectbox(f"Ambiente '{amb}' → Grupo:", [1, 2, 3], key=f"grupo_{amb}")
            ambiente_grupos[amb] = grupo

        if st.button("Ejecutar fusión y generación de centroides"):
            # --- Spatial join con tratamiento
            st.info("Spatial join de la grilla con tratamiento...")
            gdf_grid_trat = gpd.sjoin(gdf_grid, gdf_trat[[col_tratamiento, "geometry"]], how="left", predicate="intersects")
            gdf_grid_trat = gdf_grid_trat.rename(columns={col_tratamiento: "tratamiento"})

            # Dominante
            tratamiento_dominante = gdf_grid_trat.groupby(gdf_grid_trat.index)["tratamiento"].agg(
                lambda x: x.mode().iloc[0] if not x.mode().empty else np.nan
            )
            gdf_grid["tratamiento"] = gdf_grid.index.map(tratamiento_dominante)


            # --- Filtrar grilla que cayó en tratamiento
            gdf_grid_valid = gdf_grid[~gdf_grid["tratamiento"].isnull()].copy()

            # --- Spatial join con ambiente
            st.info("Spatial join de la grilla con ambiente...")
            gdf_grid_amb = gpd.sjoin(gdf_grid_valid, gdf_amb[[col_amb, "geometry"]], how="left", predicate="intersects")
            gdf_grid_amb = gdf_grid_amb.rename(columns={col_amb: "amb"})

            ambiente_dominante = gdf_grid_amb.groupby(gdf_grid_amb.index)["amb"].agg(lambda x: x.mode().iloc[0] if not x.mode().empty else np.nan)
            gdf_grid_valid["amb"] = gdf_grid_valid.index.map(ambiente_dominante)

            # --- Mapear grupo de ambiente
            gdf_grid_valid["amb_grupo"] = gdf_grid_valid["amb"].map(ambiente_grupos)

            # --- Agregar centroides
            gdf_grid_valid["xcoord"] = gdf_grid_valid.geometry.centroid.x
            gdf_grid_valid["ycoord"] = gdf_grid_valid.geometry.centroid.y

            # --- Agregar rinde
            st.info("Asignando rendimiento por nearest neighbor...")
            gdf_clean = st.session_state["gdf_clean"]
            col_rinde = st.session_state["col_rinde"]

            gdf_clean_points = gdf_clean.copy()
            gdf_clean_points["x"] = gdf_clean_points.geometry.x
            gdf_clean_points["y"] = gdf_clean_points.geometry.y

            from sklearn.neighbors import NearestNeighbors
            nn = NearestNeighbors(n_neighbors=1)
            nn.fit(gdf_clean_points[["x", "y"]].values)

            distances, indices = nn.kneighbors(gdf_grid_valid[["xcoord", "ycoord"]].values)
            gdf_grid_valid["rinde_mean"] = gdf_clean_points.iloc[indices.flatten()][col_rinde].values

            # --- Mostrar preview
            st.success("Fusión completa.")
            st.dataframe(gdf_grid_valid[["tratamiento", "amb", "amb_grupo", "rinde_mean", "xcoord", "ycoord"]].head(50))

            # --- Botón de descarga CSV
            csv_buffer = gdf_grid_valid[["tratamiento", "amb", "amb_grupo", "rinde_mean", "xcoord", "ycoord"]].to_csv(index=False).encode()
            st.download_button(
                label="📥 Descargar CSV intermedio",
                data=csv_buffer,
                file_name="grilla_intermedia.csv",
                mime="text/csv"
            )

            # Guardar en session_state para Tab 3
            st.session_state["gdf_grid_valid"] = gdf_grid_valid
# TAB 3 — Limpieza por vecinos homogéneos (nuevo enfoque sólido)

with tabs[2]:
    st.header("3️⃣ Limpieza de borde por vecinos de tratamiento distinto")

    if "gdf_grid_valid" not in st.session_state:
        st.warning("⚠️ Primero debes ejecutar la fusión en el Tab 2.")
    else:
        gdf_grid_valid = st.session_state["gdf_grid_valid"].copy()

        # --- 1️⃣ Proyectar a UTM para seguridad (no afecta geometría poligonal)
        epsg_planar = 32720
        gdf_grid_planar = gdf_grid_valid.to_crs(epsg=epsg_planar)

        # --- 2️⃣ Reset index para evitar problemas en .iloc
        gdf_grid_planar = gdf_grid_planar.reset_index(drop=True)

        # --- 3️⃣ Construir spatial index
        sindex = gdf_grid_planar.sindex

        # --- 4️⃣ Precalcular mapa de tratamiento
        mapa_tratamiento = gdf_grid_planar["tratamiento"].to_dict()

        # --- 5️⃣ Definir función de conteo de vecinos distintos
        def cuenta_vecinos_distintos(idx, geom):
            vecinos_idx_pos = list(sindex.query(geom, predicate='touches'))
            vecinos_idx_pos = [v for v in vecinos_idx_pos if v != idx]

            vecinos_distintos = 0
            for v_pos in vecinos_idx_pos:
                v_geom = gdf_grid_planar.iloc[v_pos].geometry
                if geom.touches(v_geom):
                    if mapa_tratamiento.get(v_pos, None) != mapa_tratamiento.get(idx, None):
                        vecinos_distintos += 1
            return vecinos_distintos

        # --- 6️⃣ Aplicar función
        st.info("Detectando celdas con vecinos de tratamiento distinto...")
        gdf_grid_planar["vecinos_distintos"] = gdf_grid_planar.apply(
            lambda row: cuenta_vecinos_distintos(row.name, row.geometry),
            axis=1
        )

        # --- 7️⃣ Filtrar solo celdas con 0 vecinos distintos (homogéneo)
        gdf_grid_valid_final = gdf_grid_planar[gdf_grid_planar["vecinos_distintos"] == 0].copy()

        st.success(f"Celdas finales con vecindario homogéneo: {len(gdf_grid_valid_final)}")

        # --- 8️⃣ Volver a EPSG:4326 para exportación
        gdf_grid_valid_final = gdf_grid_valid_final.to_crs(epsg=4326)
        gdf_grid_valid_final["xcoord"] = gdf_grid_valid_final.geometry.centroid.x
        gdf_grid_valid_final["ycoord"] = gdf_grid_valid_final.geometry.centroid.y

        # --- 9️⃣ Mostrar preview
        st.dataframe(gdf_grid_valid_final[["tratamiento", "amb", "amb_grupo", "rinde_mean", "xcoord", "ycoord"]].head(50))

        # --- 10️⃣ Botón de descarga CSV
        csv_buffer = gdf_grid_valid_final[["tratamiento", "amb", "amb_grupo", "rinde_mean", "xcoord", "ycoord"]].to_csv(index=False).encode()

        st.download_button(
            label="📥 Descargar CSV limpio (vecindario homogéneo)",
            data=csv_buffer,
            file_name="grilla_final_vecinos_homogeneos.csv",
            mime="text/csv"
        )
