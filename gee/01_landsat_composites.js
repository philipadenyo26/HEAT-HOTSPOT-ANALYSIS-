// =====================================================================
// Step 1. Landsat composites, NDVI and NDBI, and the scene list
// Chapter 3, Phase 2 (Data Preparation)
// Run in the GEE Code Editor. Change PROJECT to the Earth Engine project in use.
// =====================================================================
var PROJECT   = 'projects/lunar-tine-512302-q0/assets/';
var districts = ee.FeatureCollection(PROJECT + 'GAR_districts');   // 29 districts (OCHA)
var region    = districts.geometry().buffer(2000);                // 2 km buffer against edge effects
var CRS = 'EPSG:32630';                                            // WGS 84 / UTM Zone 30N
var FOLDER = 'GEE_thesis';

// Harmattan dry season: November to March
var drySeason = ee.Filter.or(ee.Filter.calendarRange(11, 12, 'month'),
                             ee.Filter.calendarRange(1, 3, 'month'));

// Cloud and shadow mask from QA_PIXEL (bit 1 dilated cloud, bit 3 cloud, bit 4 cloud shadow)
function maskL2(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0)
    .and(qa.bitwiseAnd(1 << 3).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0));
  return img.updateMask(clear);
}

// Surface reflectance scaling and the two indices reviewed in Chapter Two
function addIndices(img) {
  var sr = img.select(['SR_B4', 'SR_B5', 'SR_B6']).multiply(0.0000275).add(-0.2);
  var ndvi = sr.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');   // (NIR - Red)/(NIR + Red)
  var ndbi = sr.normalizedDifference(['SR_B6', 'SR_B5']).rename('NDBI');   // (SWIR1 - NIR)/(SWIR1 + NIR)
  return img.addBands([ndvi, ndbi]);
}

function landsatL2(start, end) {
  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
  var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');
  return l8.merge(l9)
    .filterBounds(region)
    .filterDate(start, end)
    .filter(drySeason)
    .filter(ee.Filter.lt('CLOUD_COVER', 20))
    .map(maskL2)
    .map(addIndices);
}

var p1 = landsatL2('2015-01-01', '2020-01-01');   // 2015 to 2019
var p2 = landsatL2('2020-01-01', '2024-01-01');   // 2020 to 2023 (baseline)
print('Scenes 2015-2019', p1.size(), 'Scenes 2020-2023', p2.size());

function exportImg(img, name, scale) {
  Export.image.toDrive({image: img.clip(region).toFloat(), description: name, folder: FOLDER,
    fileNamePrefix: name, region: region, scale: scale, crs: CRS, maxPixels: 1e10});
}
exportImg(p1.select('NDVI').median(), 'NDVI_2015_2019', 30);
exportImg(p1.select('NDBI').median(), 'NDBI_2015_2019', 30);
exportImg(p2.select('NDVI').median(), 'NDVI_2020_2023', 30);
exportImg(p2.select('NDBI').median(), 'NDBI_2020_2023', 30);

// Scene list for the ERA5 water vapour download (Python step)
var scenes = p1.merge(p2).map(function (img) {
  return ee.Feature(null, {
    scene_id: img.get('system:index'),
    date: img.date().format('YYYY-MM-dd'),
    hour_utc: img.date().format('HH'),
    cloud_cover: img.get('CLOUD_COVER')
  });
});
Export.table.toDrive({collection: scenes, description: 'landsat_scene_list', folder: FOLDER,
  fileFormat: 'CSV', selectors: ['scene_id', 'date', 'hour_utc', 'cloud_cover']});

Map.centerObject(region, 10);
Map.addLayer(p2.select('NDVI').median().clip(region), {min: 0, max: 0.6, palette: ['brown', 'yellow', 'green']}, 'NDVI 2020-2023');
