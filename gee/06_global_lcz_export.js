// =====================================================================
// Global LCZ map (Demuzere et al., 2022) for the study area, 100 m, UTM 30N.
// A light export (no classification): R compares it with the LCZ_100m map.
// =====================================================================
var PROJECT   = 'projects/lunar-tine-512302-q0/assets/';
var districts = ee.FeatureCollection(PROJECT + 'GAR_districts');
var region    = districts.geometry().buffer(2000);
var globalLCZ = ee.ImageCollection('RUB/RUBCLIM/LCZ/global_lcz_map/latest').mosaic().select('LCZ_Filter');
Export.image.toDrive({image: globalLCZ.clip(region).toByte(), description: 'global_LCZ_100m', folder: 'GEE_thesis',
  region: region, scale: 100, crs: 'EPSG:32630', maxPixels: 1e9});
print('Run the global_LCZ_100m task, then download global_LCZ_100m.tif into data/gee/');
