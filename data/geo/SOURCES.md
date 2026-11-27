# Map layers

`map_layers.gpkg` holds the vector layers of report Figure 2 (study areas and field sites), in EPSG:4326, cropped to the regional map and the four study-area panels. `map_places.csv` holds the towns and field bases drawn on the panels, at town centres.

| Layer | Content | Source | Licence |
|---|---|---|---|
| `countries` | Country outlines (Horn of Africa) | Natural Earth, 1:10m Admin 0 countries | Public domain |
| `admin_units` (Ethiopia) | Woreda boundaries around the two Ethiopian study areas, with Lake Tana cut out | OCHA Ethiopia administrative boundaries, admin 3 (v04, 2025), via the Humanitarian Data Exchange; Lake Tana from Natural Earth 1:10m lakes | CC BY-IGO 3.0; public domain |
| `admin_units` (Somalia) | District boundaries around the two Somali study areas | OCHA Somalia administrative boundaries, admin 2 (2023), via the Humanitarian Data Exchange | CC BY-IGO 3.0 |
| `settlements` | Settlement extents | GRID3 settlement extents: Ethiopia v3.0, Somalia v2.0 | CC BY 4.0 |
| `roads` | Trunk, primary, secondary and tertiary roads | OpenStreetMap, Humanitarian OpenStreetMap Team exports | ODbL 1.0, © OpenStreetMap contributors |
| `rivers` | Rivers in the Somali panels | OpenStreetMap waterways | ODbL 1.0, © OpenStreetMap contributors |
| `map_places.csv` | Towns and field bases | OpenStreetMap Nominatim (Ethiopia); OCHA Somalia settlements 2022 (Somalia) | ODbL 1.0; CC BY-IGO 3.0 |

The roads and rivers layers are extracts of OpenStreetMap and are made available under the Open Database Licence (ODbL) 1.0.
