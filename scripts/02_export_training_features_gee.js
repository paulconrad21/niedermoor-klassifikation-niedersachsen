// ============================================================================
// 02_export_training_features_gee.js
// Projekt Neustart: Niedermoor-Klassifikation Niedersachsen
//
// Ziel:
// Fuer jedes Trainingspolygon werden Fernerkundungs- und Kontextvariablen
// berechnet und als CSV-Tabellen nach Google Drive exportiert.
//
// Ergebnis:
// - eine Tabelle fuer Sentinel-2 Monatsfeatures
// - eine Tabelle fuer Sentinel-1 Monatsfeatures
// - eine Tabelle fuer Landsat-LST Temperaturfeatures
// - eine Tabelle fuer JRC Wasser-/Oberflaechenwasserfeatures
//
// Wichtig:
// Dieses Skript exportiert nur Werte fuer Trainingspolygone.
// Die flaechendeckende Karte fuer ganz Niedersachsen kommt erst spaeter,
// nachdem klar ist, welche Features wirklich gebraucht werden.
// ============================================================================


// ============================================================================
// 1. Einstellungen
// ============================================================================

// TODO: Nach dem Upload in Google Earth Engine hier den Asset-Pfad anpassen.
// Beispiel:
// var trainingAsset = 'projects/geo-erfassung-data/assets/training_polygons_clean';
var trainingAsset = 'projects/geo-erfassung-data/assets/training_polygons_clean';

// Untersuchungsjahr.
var year = 2025;

// Zielordner in Google Drive.
// Der Ordner wird von GEE beim Export angelegt, falls er noch nicht existiert.
var driveFolder = 'GEE_exports_neustart_2025';

// Skalen:
// Sentinel-2 hat 10-m- und 20-m-Baender. Da wir auch Red-Edge/SWIR nutzen,
// ist 20 m fuer gemischte S2-Features die konservativere Wahl.
var s2Scale = 20;

// Sentinel-1 liegt in GEE ueblicherweise mit 10 m Pixelgroesse vor.
var s1Scale = 10;

// Landsat LST und JRC Surface Water sind 30-m-Produkte.
var landsatScale = 30;
var jrcScale = 30;

// Bei grossen oder komplexen Polygonen kann tileScale helfen,
// Speicherprobleme in GEE zu reduzieren.
var tileScale = 4;

// Diese Eigenschaften sollen in allen Exporttabellen erhalten bleiben.
// Sie kommen aus Skript 01 bzw. dem GEE-Shapefile.
var idProperties = ['poly_id', 'is_moor', 'cls_bin', 'cls_det'];


// ============================================================================
// 2. Trainingspolygone laden
// ============================================================================

var training = ee.FeatureCollection(trainingAsset);
var region = training.geometry();

print('Trainingspolygone:', training.size());
print('Beispielpolygon:', training.first());

Map.centerObject(training, 8);
Map.addLayer(training, {color: 'orange'}, 'Trainingspolygone', true);


// ============================================================================
// 3. Kleine Hilfsfunktionen
// ============================================================================

// Erzeugt eine Monatsliste von 1 bis 12.
var months = ee.List.sequence(1, 12);

// Erzeugt eine Monatskennung wie "m01", "m02", ...
function monthLabel(month) {
  month = ee.Number(month);
  return ee.String('m').cat(month.format('%02d'));
}

// Benennt alle Baender eines Bildes mit einem Prefix um.
// Beispiel: NDVI -> S2_m03_NDVI
function prefixBands(image, prefix) {
  var oldNames = image.bandNames();
  var newNames = oldNames.map(function(name) {
    return ee.String(prefix).cat('_').cat(ee.String(name));
  });

  return image.rename(newNames);
}

// Leeres, komplett maskiertes Bild mit festen Bandnamen.
// Das ist wichtig fuer Monate ohne gueltige Szenen.
// Dadurch bleibt die Tabellenstruktur ueber alle Monate stabil.
function emptyImage(bandNames) {
  var zeros = ee.Image.constant(ee.List.repeat(0, bandNames.length()));
  return zeros.rename(bandNames).updateMask(ee.Image(0));
}

// Exportiert ein Bild als breite Tabelle: eine Zeile pro Polygon,
// viele Spalten fuer die Featurewerte.
function exportImageByPolygons(image, reducer, scale, description) {
  var table = image.reduceRegions({
    collection: training,
    reducer: reducer,
    scale: scale,
    crs: 'EPSG:25832',
    tileScale: tileScale
  });

  // Geometrie wird fuer CSV nicht benoetigt.
  // poly_id verbindet die Tabelle spaeter wieder mit den Polygonen in R.
  table = table.map(function(feature) {
    return feature.setGeometry(null);
  });

  Export.table.toDrive({
    collection: table,
    description: description,
    folder: driveFolder,
    fileNamePrefix: description,
    fileFormat: 'CSV'
  });
}


// ============================================================================
// 4. Sentinel-2 vorbereiten
// ============================================================================

// Sentinel-2 Surface Reflectance.
var s2Raw = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(region)
  .filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year + 1, 1, 1))
  // Nicht zu streng filtern, weil wir spaeter mit SCL maskieren.
  // Ein strenger Filter kann in wolkigen Monaten zu Datenluecken fuehren.
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 85));

print('Sentinel-2 Szenen im Jahr:', s2Raw.size());

// Wolkenmaske ueber die Scene Classification Layer (SCL).
// Entfernt: Schatten, Wolken, Cirrus, Schnee/Eis und fehlerhafte Pixel.
function maskS2Clouds(image) {
  var scl = image.select('SCL');

  var clearMask = scl.neq(1)   // saturated/defective
    .and(scl.neq(3))           // cloud shadow
    .and(scl.neq(8))           // cloud medium probability
    .and(scl.neq(9))           // cloud high probability
    .and(scl.neq(10))          // thin cirrus
    .and(scl.neq(11));         // snow/ice

  return image.updateMask(clearMask);
}

// Rohbaender umbenennen und wichtige Indizes berechnen.
function addS2Features(image) {
  image = maskS2Clouds(image);

  // S2-SR ist mit Faktor 10000 skaliert.
  var scaled = image.select(
    ['B2',   'B3',    'B4',  'B5',  'B6',  'B7',  'B8', 'B8A',   'B11',   'B12'],
    ['blue', 'green', 'red', 're1', 're2', 're3', 'nir','nir8a', 'swir1', 'swir2']
  ).multiply(0.0001);

  var ndvi = scaled.normalizedDifference(['nir', 'red']).rename('NDVI');
  var ndmi = scaled.normalizedDifference(['nir', 'swir1']).rename('NDMI');
  var ndwi = scaled.normalizedDifference(['green', 'nir']).rename('NDWI');
  var mndwi = scaled.normalizedDifference(['green', 'swir1']).rename('MNDWI');
  var ndre = scaled.normalizedDifference(['nir8a', 're1']).rename('NDRE');
  var nbr = scaled.normalizedDifference(['nir', 'swir2']).rename('NBR');
  var nbr2 = scaled.normalizedDifference(['swir1', 'swir2']).rename('NBR2');

  // EVI ist oft robuster als NDVI bei dichter Vegetation.
  var evi = scaled.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
    {
      NIR: scaled.select('nir'),
      RED: scaled.select('red'),
      BLUE: scaled.select('blue')
    }
  ).rename('EVI');

  return scaled
    .addBands([ndvi, ndmi, ndwi, mndwi, ndre, nbr, nbr2, evi])
    .copyProperties(image, image.propertyNames());
}

var s2Features = s2Raw.map(addS2Features);

var s2BandNames = [
  'blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir8a', 'swir1', 'swir2',
  'NDVI', 'NDMI', 'NDWI', 'MNDWI', 'NDRE', 'NBR', 'NBR2', 'EVI'
];

// Monatskomposit fuer Sentinel-2.
function s2MonthlyComposite(month) {
  month = ee.Number(month);

  var start = ee.Date.fromYMD(year, month, 1);
  var end = start.advance(1, 'month');
  var monthlyCollection = s2Features.filterDate(start, end);

  var image = ee.Image(
    ee.Algorithms.If(
      monthlyCollection.size().gt(0),
      monthlyCollection.select(s2BandNames).median(),
      emptyImage(s2BandNames)
    )
  );

  return prefixBands(image, ee.String('S2_').cat(monthLabel(month)));
}

var s2MonthlyStack = ee.ImageCollection.fromImages(
  months.map(s2MonthlyComposite)
).toBands();

// toBands erzeugt zusaetzliche Index-Prefixe. Deshalb werden die Bandnamen bereinigt.
s2MonthlyStack = s2MonthlyStack.rename(
  s2MonthlyStack.bandNames().map(function(name) {
    name = ee.String(name);
    return name.replace('^[0-9]+_', '');
  })
);

// Nicht mit print() oder Map.addLayer() auf den kompletten S2-Stack zugreifen.
// Das kann in GEE schon vor dem Export ein Speicherlimit ausloesen.
// Erwartung: 12 Monate * 18 S2-Features = 216 Baender.
print('S2 Feature-Baender erwartet:', 216);


// ============================================================================
// 5. Sentinel-1 vorbereiten
// ============================================================================

var s1Raw = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(region)
  .filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year + 1, 1, 1))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));

print('Sentinel-1 Szenen im Jahr:', s1Raw.size());

// Sentinel-1 wird in GEE als dB-Wert bereitgestellt.
// Wir behalten ASC und DESC getrennt, weil unterschiedliche Blickrichtungen
// nicht einfach als identisch behandelt werden sollten.
function addS1Features(image) {
  var vv = image.select('VV').rename('VV');
  var vh = image.select('VH').rename('VH');
  var vvMinusVh = vv.subtract(vh).rename('VV_minus_VH');

  return ee.Image.cat([vv, vh, vvMinusVh])
    .copyProperties(image, image.propertyNames());
}

var s1Features = s1Raw.map(addS1Features);
var s1BandNames = ['VV', 'VH', 'VV_minus_VH'];

function s1MonthlyComposite(month, orbit) {
  month = ee.Number(month);

  var start = ee.Date.fromYMD(year, month, 1);
  var end = start.advance(1, 'month');

  var monthlyCollection = s1Features
    .filterDate(start, end)
    .filter(ee.Filter.eq('orbitProperties_pass', orbit));

  var image = ee.Image(
    ee.Algorithms.If(
      monthlyCollection.size().gt(0),
      monthlyCollection.select(s1BandNames).median(),
      emptyImage(s1BandNames)
    )
  );

  var prefix = ee.String('S1_')
    .cat(monthLabel(month))
    .cat('_')
    .cat(orbit);

  return prefixBands(image, prefix);
}

var s1AscStack = ee.ImageCollection.fromImages(
  months.map(function(month) {
    return s1MonthlyComposite(month, 'ASCENDING');
  })
).toBands();

var s1DescStack = ee.ImageCollection.fromImages(
  months.map(function(month) {
    return s1MonthlyComposite(month, 'DESCENDING');
  })
).toBands();

var s1MonthlyStack = s1AscStack.addBands(s1DescStack);

s1MonthlyStack = s1MonthlyStack.rename(
  s1MonthlyStack.bandNames().map(function(name) {
    name = ee.String(name);
    return name
      .replace('^[0-9]+_', '')
      .replace('ASCENDING', 'ASC')
      .replace('DESCENDING', 'DESC');
  })
);

print('S1 Feature-Baender:', s1MonthlyStack.bandNames().size());


// ============================================================================
// 6. Landsat Land Surface Temperature vorbereiten
// ============================================================================

// Landsat 8 und 9 Level 2 enthalten Surface Temperature in ST_B10.
var landsat8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
var landsat9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');

var landsatRaw = landsat8.merge(landsat9)
  .filterBounds(region)
  .filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year + 1, 1, 1))
  // L2SP bedeutet: Surface Temperature ist vorhanden.
  .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'));

print('Landsat L2SP Szenen im Jahr:', landsatRaw.size());

function maskLandsatL2(image) {
  var qa = image.select('QA_PIXEL');

  // QA_PIXEL Bits:
  // 0 Fill, 1 Dilated Cloud, 2 Cirrus, 3 Cloud,
  // 4 Cloud Shadow, 5 Snow
  var mask = qa.bitwiseAnd(1 << 0).eq(0)
    .and(qa.bitwiseAnd(1 << 1).eq(0))
    .and(qa.bitwiseAnd(1 << 2).eq(0))
    .and(qa.bitwiseAnd(1 << 3).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0))
    .and(qa.bitwiseAnd(1 << 5).eq(0));

  return image.updateMask(mask);
}

function addLandsatLST(image) {
  image = maskLandsatL2(image);

  // Skalierung laut Landsat Collection 2 Level 2:
  // Kelvin = DN * 0.00341802 + 149.0
  // Celsius = Kelvin - 273.15
  var lstC = image.select('ST_B10')
    .multiply(0.00341802)
    .add(149.0)
    .subtract(273.15)
    .rename('LST_C');

  return lstC.copyProperties(image, image.propertyNames());
}

var landsatLST = landsatRaw.map(addLandsatLST);

function seasonalLST(startMonth, endMonth, prefix) {
  var start = ee.Date.fromYMD(year, startMonth, 1);
  var end = ee.Date.fromYMD(year, endMonth, 1).advance(1, 'month');
  var collection = landsatLST.filterDate(start, end);

  var mean = collection.mean().rename(prefix + '_mean_C');
  var min = collection.min().rename(prefix + '_min_C');
  var max = collection.max().rename(prefix + '_max_C');
  var amp = max.subtract(min).rename(prefix + '_amplitude_C');

  return ee.Image.cat([mean, min, max, amp]);
}

var landsatLSTStack = ee.Image.cat([
  seasonalLST(3, 5, 'LST_spring'),
  seasonalLST(6, 8, 'LST_summer'),
  seasonalLST(3, 9, 'LST_growing')
]);

print('Landsat LST Feature-Baender:', landsatLSTStack.bandNames());


// ============================================================================
// 7. JRC Global Surface Water vorbereiten
// ============================================================================

// Dieser Datensatz beschreibt Oberflaechenwasser-Dynamik.
// Er ist kein direkter Grundwasser- oder Moorboden-Datensatz,
// aber ein sinnvoller Wasser-Kontext.
var jrc = ee.Image('JRC/GSW1_4/GlobalSurfaceWater');

var jrcWater = jrc.select(
  ['occurrence', 'seasonality', 'recurrence', 'max_extent'],
  ['JRC_occurrence', 'JRC_seasonality', 'JRC_recurrence', 'JRC_max_extent']
);

// Distanz zum historisch maximal beobachteten Oberflaechenwasser.
// Radius 10 km: groesser ist oft teurer und fuer Trainingsfeatures selten noetig.
var maxWaterMask = jrc.select('max_extent').eq(1).selfMask();
// GEE erlaubt nur Distanz-Kernel bis 512 Pixel.
// 5000 m ist hier bewusst konservativ gewaehlt, damit der Export stabil laeuft.
var maxWaterDistanceM = 5000;

var distanceToWater = maxWaterMask
  .distance(ee.Kernel.euclidean(maxWaterDistanceM, 'meters'))
  .unmask(maxWaterDistanceM)
  .rename('JRC_dist_max_extent_m');

var jrcStack = jrcWater.addBands(distanceToWater);

print('JRC Feature-Baender:', jrcStack.bandNames());

Map.addLayer(
  jrc.select('occurrence'),
  {min: 0, max: 100, palette: ['white', 'lightblue', 'blue']},
  'Kontrolle: JRC Wasser occurrence',
  false
);


// ============================================================================
// 8. Tabellenexporte starten
// ============================================================================

// Sentinel-2:
// Median pro Polygon, weil einzelne Ausreisser weniger stark wirken.
exportImageByPolygons(
  s2MonthlyStack,
  ee.Reducer.median(),
  s2Scale,
  'training_features_s2_2025'
);

// Sentinel-1:
// Ebenfalls Median pro Polygon.
exportImageByPolygons(
  s1MonthlyStack,
  ee.Reducer.median(),
  s1Scale,
  'training_features_s1_2025'
);

// Landsat LST:
// Mittelwert pro Polygon ist fuer saisonale Temperaturfeatures gut interpretierbar.
exportImageByPolygons(
  landsatLSTStack,
  ee.Reducer.mean(),
  landsatScale,
  'training_features_landsat_lst_2025'
);

// JRC Wasser:
// Mittelwert pro Polygon fuer occurrence/seasonality/recurrence.
// Fuer Distanz ergibt der Mittelwert eine ungefaehre mittlere Wassernaehe.
exportImageByPolygons(
  jrcStack,
  ee.Reducer.mean(),
  jrcScale,
  'training_features_jrc_water'
);

print('Exports wurden vorbereitet. Bitte rechts im Tasks-Tab alle vier Exporte starten.');
