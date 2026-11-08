// =====================================================================
// WorldPop 2020 (100 m, unconstrained) for Ghana, exported from GEE
// Replaces the manual WorldPop download. Used in Step 6 (heat risk).
// Exports: total population and the vulnerable age group (under 10 and 65+),
// both as people per cell, in WorldPop's native grid.
// =====================================================================
var PROJECT   = 'projects/lunar-tine-512302-q0/assets/';
var districts = ee.FeatureCollection(PROJECT + 'GAR_districts');
var region    = districts.geometry().buffer(5000);
var FOLDER = 'GEE_thesis';

var ageSex = ee.ImageCollection('WorldPop/GP/100m/pop_age_sex')
  .filter(ee.Filter.eq('country', 'GHA')).filter(ee.Filter.eq('year', 2020)).first();

var total = ageSex.select('population').rename('pop_total');
var vulBands = ['M_0', 'M_1', 'M_5', 'M_65', 'M_70', 'M_75', 'M_80',
                'F_0', 'F_1', 'F_5', 'F_65', 'F_70', 'F_75', 'F_80'];
var vulnerable = ageSex.select(vulBands).reduce(ee.Reducer.sum()).rename('pop_vulnerable');

print('Total population in region', total.reduceRegion({reducer: ee.Reducer.sum(),
  geometry: districts.geometry(), scale: 100, maxPixels: 1e10}));

var p = ageSex.select('population').projection().getInfo();   // client-side CRS and transform
Export.image.toDrive({image: total.clip(region).toFloat(), description: 'worldpop_total_2020',
  folder: FOLDER, region: region, crs: p.crs, crsTransform: p.transform, maxPixels: 1e10});
Export.image.toDrive({image: vulnerable.clip(region).toFloat(), description: 'worldpop_vulnerable_2020',
  folder: FOLDER, region: region, crs: p.crs, crsTransform: p.transform, maxPixels: 1e10});
