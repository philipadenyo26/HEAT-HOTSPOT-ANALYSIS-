// =====================================================================
// Step 3. LST retrieval with the Single-Channel algorithm (Band 10)
// Chapter 3, Phase 3 (Objective 1)
// Inputs: Level-1 TOA brightness temperature (Band 10), Level-2 surface
// reflectance for NDVI and the QA mask, ERA5 water vapour per scene.
// =====================================================================
var PROJECT   = 'projects/lunar-tine-512302-q0/assets/';
var districts = ee.FeatureCollection(PROJECT + 'GAR_districts');
var wTable    = ee.FeatureCollection(PROJECT + 'era5_tcwv_scenes');   // columns: scene_id, date, w (g/cm2)
var region    = districts.geometry().buffer(2000);
var CRS = 'EPSG:32630';
var FOLDER = 'GEE_thesis';
var drySeason = ee.Filter.or(ee.Filter.calendarRange(11, 12, 'month'),
                             ee.Filter.calendarRange(1, 3, 'month'));
var wMean = ee.Number(wTable.aggregate_mean('w'));

// Water vapour for the scene date (falls back to the mean if a date is missing)
function getW(img) {
  var d = img.date().format('YYYY-MM-dd');
  var f = wTable.filter(ee.Filter.eq('date', d));
  return ee.Number(ee.Algorithms.If(f.size().gt(0), f.first().get('w'), wMean));
}

// Pair each Level-2 scene with its Level-1 TOA scene (same system:index), per sensor
function paired(l2id, toaid, start, end) {
  var filt = function (c) {
    return ee.ImageCollection(c).filterBounds(region).filterDate(start, end)
      .filter(drySeason).filter(ee.Filter.lt('CLOUD_COVER', 20));
  };
  var join = ee.Join.inner().apply(filt(l2id), filt(toaid),
    ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'}));
  return ee.ImageCollection(join.map(function (p) {
    return singleChannel(ee.Image(p.get('primary')), ee.Image(p.get('secondary')));
  }));
}

function singleChannel(l2, toa) {
  // Cloud and shadow mask (QA_PIXEL bits 1, 3, 4)
  var qa = l2.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0).and(qa.bitwiseAnd(1 << 3).eq(0)).and(qa.bitwiseAnd(1 << 4).eq(0));

  // NDVI, vegetation proportion and emissivity (NDVI Threshold Method, Section 3.7)
  var sr   = l2.select(['SR_B4', 'SR_B5']).multiply(0.0000275).add(-0.2);
  var ndvi = sr.normalizedDifference(['SR_B5', 'SR_B4']);
  var pv   = ndvi.subtract(0.2).divide(0.3).clamp(0, 1).pow(2).rename('PV');           // Equation 3.2
  var mixed = pv.multiply(0.99).add(ee.Image(1).subtract(pv).multiply(0.97));
  var eps  = ee.Image(0.97).where(ndvi.gt(0.5), 0.99)
               .where(ndvi.gte(0.2).and(ndvi.lte(0.5)), mixed).rename('EMIS');

  // Band 10 brightness temperature (K) and at-sensor radiance by inverting Planck's law
  var T  = toa.select('B10');
  var K1 = ee.Number(toa.get('K1_CONSTANT_BAND_10'));
  var K2 = ee.Number(toa.get('K2_CONSTANT_BAND_10'));
  var L  = T.expression('K1 / (exp(K2 / T) - 1)', {K1: K1, K2: K2, T: T});

  // Atmospheric functions from water vapour w (Jimenez-Munoz et al., 2014, Landsat 8 Band 10)
  var w = getW(l2);
  var psi1 = w.pow(2).multiply(0.04019).add(w.multiply(0.02916)).add(1.01523);
  var psi2 = w.pow(2).multiply(-0.38333).subtract(w.multiply(1.50294)).add(0.20324);
  var psi3 = w.pow(2).multiply(0.00918).add(w.multiply(1.36072)).subtract(0.27514);

  var bGamma = 1324;                                   // K, Band 10
  var gamma = T.pow(2).divide(L.multiply(bGamma));
  var delta = T.subtract(T.pow(2).divide(bGamma));
  var lst = gamma.multiply(L.multiply(psi1).add(psi2).divide(eps).add(psi3))
              .add(delta).subtract(273.15).rename('LST');   // degrees Celsius

  return lst.addBands([pv, eps]).updateMask(clear)
            .set('system:time_start', l2.get('system:time_start'))
            .set('w', w);
}

function lstPeriod(start, end) {
  return paired('LANDSAT/LC08/C02/T1_L2', 'LANDSAT/LC08/C02/T1_TOA', start, end)
    .merge(paired('LANDSAT/LC09/C02/T1_L2', 'LANDSAT/LC09/C02/T1_TOA', start, end));
}
var p1 = lstPeriod('2015-01-01', '2020-01-01');
var p2 = lstPeriod('2020-01-01', '2024-01-01');

function exportImg(img, name) {
  Export.image.toDrive({image: img.clip(region).toFloat(), description: name, folder: FOLDER,
    fileNamePrefix: name, region: region, scale: 30, crs: CRS, maxPixels: 1e10});
}
exportImg(p1.select('LST').median(), 'LST_2015_2019');
exportImg(p2.select('LST').median(), 'LST_2020_2023');
exportImg(p2.select('PV').median(), 'PV_2020_2023');

// Per-scene median LST over the region, for the temporal check in the EDA
var perScene = p1.merge(p2).map(function (img) {
  var m = img.select('LST').reduceRegion({reducer: ee.Reducer.median(), geometry: region,
    scale: 90, maxPixels: 1e9, tileScale: 4});
  return ee.Feature(null, {date: img.date().format('YYYY-MM-dd'), w: img.get('w'), lst_median: m.get('LST')});
});
Export.table.toDrive({collection: perScene, description: 'lst_per_scene', folder: FOLDER,
  fileFormat: 'CSV', selectors: ['date', 'w', 'lst_median']});

Map.centerObject(region, 10);
Map.addLayer(p2.select('LST').median().clip(region), {min: 28, max: 48,
  palette: ['313695', '74add1', 'fed976', 'feb24c', 'f03b20', 'bd0026']}, 'LST 2020-2023');
