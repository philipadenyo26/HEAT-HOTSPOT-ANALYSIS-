// =====================================================================
// Independent LST check: USGS Landsat Collection 2 Level-2 surface temperature (ST_B10),
// dry-season median 2020-2023, same scenes, masks and grid as LST_2020_2023.
// A light export (no Single-Channel computation). R compares it with the Single-Channel LST.
// =====================================================================
var PROJECT   = 'projects/lunar-tine-512302-q0/assets/';
var districts = ee.FeatureCollection(PROJECT + 'GAR_districts');
var region    = districts.geometry().buffer(2000);
var dry = ee.Filter.or(ee.Filter.calendarRange(11, 12, 'month'), ee.Filter.calendarRange(1, 3, 'month'));
function st(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0).and(qa.bitwiseAnd(1 << 3).eq(0)).and(qa.bitwiseAnd(1 << 4).eq(0));
  return img.select('ST_B10').multiply(0.00341802).add(149.0).subtract(273.15)   // USGS scale and offset, to deg C
            .rename('ST_USGS').updateMask(clear);
}
var col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
  .filterBounds(region).filterDate('2020-01-01', '2024-01-01').filter(dry)
  .filter(ee.Filter.lt('CLOUD_COVER', 20)).map(st);
Export.image.toDrive({image: col.median().clip(region).toFloat(), description: 'USGS_ST_2020_2023',
  folder: 'GEE_thesis', region: region, scale: 30, crs: 'EPSG:32630', maxPixels: 1e10});
print('Run the USGS_ST_2020_2023 task, then download USGS_ST_2020_2023.tif into data/gee/');
