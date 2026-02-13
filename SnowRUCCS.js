/**
 * SnowRUCCS - Snow Recognition Under Cloud Cover System
 * 
 * An algorithm for reconstructing snow information under clouds and cloud shadows
 * using Sentinel-2 imagery based on FMask 4.0 cloud detection.
 * 
 * Features:
 *   - Two-pass cloud detection algorithm (PCP + probability-based refinement)
 *   - Improved cloud shadow detection using FMask fillMinima and terrain shadow masking
 *   - Elevation-stratified snow reconstruction under cloud-covered areas
 *   - Multi-threshold snow detection with NDSI, NIR, and spectral ratios
 * 
 * Data Sources:
 *   - Sentinel-2 Level-1C (Cirrus band B10)
 *   - Sentinel-2 Level-2A (Surface Reflectance)
 *   - JRC Global Surface Water v1.3
 *   - SRTM DEM 30m
 *   - ESA WorldCover v200
 * 
 * Usage:
 *   Run this script in Google Earth Engine Code Editor.
 *   Define your study area as 'roi' before execution.
 * 
 * Author: SnowRUCCS Contributors
 * License: MIT
 */

// ============================================================
// 0. CONFIGURATION - Modify these parameters for your study
// ============================================================

// Date range for image acquisition
var s_date = ee.Date('2020-03-21');
var e_date = ee.Date('2020-03-22');

// Study area - Replace with your own FeatureCollection or Geometry
// Option 1: Use a FeatureCollection asset
// var roi = ee.FeatureCollection("users/your_username/your_study_area").geometry();

// Option 2: Define a polygon manually
// var roi = ee.Geometry.Polygon([
//   [[lon1, lat1], [lon2, lat2], [lon3, lat3], [lon4, lat4], [lon1, lat1]]
// ]);

// For demonstration, using a sample study area
var table = ee.FeatureCollection("users/ycqysl/BBH");
var roi = table.geometry();
Map.centerObject(roi, 10);
Map.addLayer(roi, {}, 'Study Area');

// ============================================================
// 1. DATA LOADING AND PREPROCESSING
// ============================================================

// Load Sentinel-2 Level-1C for Cirrus band
var S21C = ee.ImageCollection("COPERNICUS/S2")
    .filterDate(s_date, e_date)
    .filterBounds(roi)
    .select("B10");

// Load Sentinel-2 Level-2A for surface reflectance bands
var S22A = ee.ImageCollection("COPERNICUS/S2_SR")
    .filterDate(s_date, e_date)
    .filterBounds(roi)
    .select("B1", "B2", "B3", "B4", "B8", "B9", "B11", "B12");

// Join L1C and L2A products by acquisition time
var innerJoin = ee.Join.inner();
var filterTimeEq = ee.Filter.equals({
  leftField: 'system:time_start',
  rightField: 'system:time_start'
});

var innerJoinedSentinel2 = innerJoin.apply(S21C, S22A, filterTimeEq);
print('Joined Sentinel-2 products:', innerJoinedSentinel2);

// Merge bands from both products
var joinedSentinel2 = innerJoinedSentinel2.map(function(feature) {
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'));
});
joinedSentinel2 = ee.ImageCollection(joinedSentinel2);

// Load ancillary datasets
var GSWO = ee.Image("JRC/GSW1_3/GlobalSurfaceWater");
var SRTM = ee.Image("USGS/SRTMGL1_003");
var worldcover = ee.ImageCollection("ESA/WorldCover/v200").first();
var worldWideGeom = ee.Geometry.Rectangle([-180, -90, 180, 90], 'EPSG:4326', false);

// Visualization parameters
var VisParam = {
  opacity: 1,
  bands: ["NIR", "Red", "Green"],
  min: 0.0674,
  max: 0.8187,
  gamma: 1.364
};

/**
 * Preprocess Sentinel-2 imagery
 * - Rename bands to standard names
 * - Scale reflectance values
 * - Add terrain and land cover auxiliary bands
 */
function RenameAddBandsS2(img) {
  // Auxiliary data layers
  var water = GSWO.select('occurrence');
  var slope = ee.Terrain.slope(SRTM).rename('slope');
  var aspect = ee.Terrain.aspect(SRTM).rename('aspect');
  var elevation = SRTM.select('elevation').rename('elevation');
  var bare = worldcover.select('Map').rename('bare');

  // Rename spectral bands
  var renamed = img.select(
    ['B2', 'B3', 'B4', 'B8', 'B10', 'B11', 'B12'],
    ['Blue', 'Green', 'Red', 'NIR', 'Cirrus', 'SWIR1', 'SWIR2']
  );
  
  // Scale to reflectance (0-1)
  var scaled = ee.Image(renamed).divide(10000);
  
  // Add all auxiliary bands
  var addBands = scaled
    .addBands(water)
    .addBands(slope)
    .addBands(aspect)
    .addBands(elevation)
    .addBands(bare);

  // Preserve metadata
  var time_start = img.get("system:time_start");
  var index = img.get('system:index');
  var out = addBands
    .set("system:time_start", time_start)
    .set('system:index', index);
  
  return out.copyProperties(img);
}

// Apply preprocessing
var S2 = joinedSentinel2.map(RenameAddBandsS2).sort("system:time_start");
print('Preprocessed Sentinel-2:', S2);
Map.addLayer(S2, VisParam, 'Sentinel-2 RGB');

// Mosaic multiple tiles if needed
var S2mosaic = S2.mosaic().clip(roi);

// ============================================================
// 2. SPECTRAL INDICES
// ============================================================

/** Normalized Difference Snow Index (NDSI) */
function NDSI(img) {
  return img.normalizedDifference(['Green', 'SWIR1']);
}

/** Normalized Difference Vegetation Index (NDVI) */
function NDVI(img) {
  return img.normalizedDifference(['NIR', 'Red']);
}

/** Normalized Difference Built-up Index (NDBI) */
function NDBI(img) {
  return img.normalizedDifference(['SWIR1', 'NIR']);
}

// ============================================================
// 3. PASS ONE - Potential Cloud Pixel (PCP) Detection
// ============================================================

/**
 * Basic spectral tests for cloud detection
 * Based on FMask algorithm criteria
 */
function BasictestS2(img) {
  var basic = img.select('SWIR2').gt(0.03)
    .and(NDSI(img).lt(0.8))
    .and(NDVI(img).lt(0.8));
  return basic;
}

/**
 * Whiteness test
 * Clouds appear white (equal reflectance across visible bands)
 */
function Whitetest(img) {
  var ADD = img.select('Blue').add(img.select('Green')).add(img.select('Red'));
  var mean = ADD.reduce(ee.Reducer.mean());
  var whiteness = img.expression(
    '(abs(B1 - mean) + abs(B2 - mean) + abs(B3 - mean)) / mean', {
      'B1': img.select('Blue'),
      'B2': img.select('Green'),
      'B3': img.select('Red'),
      'mean': mean
    });
  return whiteness.lt(0.7).multiply(whiteness).rename('whiteness');
}

/**
 * Haze Optimized Transformation (HOT) test
 * Detects haze and thin clouds
 */
function HOT(img) {
  var hot = img.expression('B1 - 0.5 * B3 - 0.08', {
    'B1': img.select('Blue'),
    'B3': img.select('Red')
  });
  return hot.gt(0).rename('HOT');
}

/**
 * NIR/SWIR1 ratio test
 * Distinguishes clouds from snow
 */
function B4B5(img) {
  var ratio = img.select('NIR').divide(img.select('SWIR1'));
  return ratio.gt(0.75).rename('Ratio4_5');
}

/**
 * Cirrus band probability
 */
function Cirrus_p(img) {
  var p = img.select('Cirrus').divide(4);
  return p.where(p.lt(0), 0);
}

/**
 * Potential Cloud Pixel (PCP) - Combined test
 */
function PCPS2(img) {
  var pcp = BasictestS2(img).eq(1)
    .and(Whitetest(img).lt(0.7))
    .and(HOT(img).eq(1))
    .and(B4B5(img).eq(1))
    .and(Cirrus_p(img).gt(0.01));
  return pcp.rename('PCP');
}

var pcpcover = S2.map(PCPS2).first();
Map.addLayer(pcpcover, {palette: ['white', 'yellow']}, 'PCP Mask');

// ============================================================
// 4. WATER AND LAND MASKS
// ============================================================

/**
 * Water body detection
 * Combines spectral tests with JRC water occurrence data
 */
function Water(img) {
  var condi1 = NDVI(img).lt(0.01).and(img.select('NIR').lt(0.11));
  var condi2 = NDVI(img).lt(0.1).and(NDVI(img).gt(0)).and(img.select('NIR').lt(0.05));
  var water0 = condi1.eq(1).or(condi2.eq(1));
  var wat = water0.eq(1);
  var occur_gs = img.select('occurrence').unmask().neq(0);
  return wat.and(occur_gs.eq(1)).rename('water');
}

/**
 * Land mask (non-water, non-cloud pixels)
 */
function LandS2(img) {
  return Water(img).eq(0).and(PCPS2(img).eq(0)).rename('land');
}

// ============================================================
// 5. PASS TWO - Cloud Probability Refinement
// ============================================================

/**
 * HOT-based cloud probability for land pixels
 */
function lHOTS2(img) {
  var HOTimg = img.expression('B1 - 0.5 * B3 - 0.08', {
    'B1': img.select('Blue'),
    'B3': img.select('Red')
  }).rename('HOT');
  
  var clear_HOT = HOTimg;
  var reducer1 = ee.Reducer.percentile([17.5]);
  var reducer2 = ee.Reducer.percentile([82.5]);
  
  var hot_l = ee.Number(clear_HOT.reduceRegion({
    reducer: reducer1,
    geometry: roi,
    maxPixels: 10e13
  }).get('HOT')).subtract(0.04);
  
  var hot_h = ee.Number(clear_HOT.reduceRegion({
    reducer: reducer2,
    geometry: roi,
    maxPixels: 10e13
  }).get('HOT')).add(0.04);
  
  var lb_p = clear_HOT.expression('(hot - hotl) / (hoth - hotl)', {
    hot: clear_HOT,
    hotl: hot_l,
    hoth: hot_h
  });
  
  var result = lb_p.where(lb_p.gt(1), 1);
  return result.where(result.lt(0), 0).rename('lHOT');
}

/**
 * Spectral variability probability
 * Low variability indicates potential cloud pixels
 */
function Vari_pS2(img) {
  // Handle saturated green band for NDSI
  var gre_satu = img.select('Green').updateMask(LandS2(img)).gte(1);
  var ndsi = NDSI(img).updateMask(LandS2(img));
  var condition1 = gre_satu.and(ndsi.lt(0)).eq(1);
  var modify_ndsi = ndsi.where(condition1, 0).abs();

  // Handle saturated red band for NDVI
  var red_satu = img.select('Red').updateMask(LandS2(img)).gte(1);
  var ndvi = NDVI(img).updateMask(LandS2(img));
  var condition2 = red_satu.and(ndvi.gt(0)).eq(1);
  var modify_ndvi = ndvi.where(condition2, 0).abs();

  var ndbi = NDBI(img).updateMask(LandS2(img));
  var maxi = modify_ndsi.max(modify_ndvi).max(ndbi).max(Whitetest(img));

  var lvari_P = maxi.expression('1 - max', {max: maxi});
  var result = lvari_P.where(lvari_P.gt(1), 1);
  return result.where(result.lt(0), 0).rename('Vari_pS2');
}

/**
 * Combined cloud probability
 */
function S2Cloud_p(img) {
  var prob = Vari_pS2(img).multiply(lHOTS2(img))
    .add(img.select('Cirrus').multiply(0.5))
    .add(0.2);
  return prob.rename('Cloud_p');
}

var S2_lcp = S2.map(S2Cloud_p).first();
Map.addLayer(S2_lcp, {min: 0, max: 1, palette: ['white', 'gray', 'black']}, 'Cloud Probability');

// ============================================================
// 6. FINAL CLOUD MASK
// ============================================================

/**
 * Generate final cloud mask with morphological operations
 */
function PotentialCloud_S2(img) {
  var cloud1 = PCPS2(img).eq(1).and(Water(img).eq(0));
  var cloud2 = cloud1.eq(1).or(S2Cloud_p(img).gt(0.9)).unmask(1);
  
  // Morphological opening to remove noise
  var open = cloud2
    .focal_min({radius: 1, kernelType: 'circle', units: 'pixels', iterations: 1})
    .focal_max({radius: 1, kernelType: 'circle', units: 'pixels', iterations: 1});
  
  return open.rename('Cloud');
}

var TCS2 = S2.map(PotentialCloud_S2).first();
print('Final Cloud Mask:', TCS2);
Map.addLayer(TCS2, {palette: ['black', 'yellow']}, 'Cloud Mask');

// ============================================================
// 7. CLOUD SHADOW DETECTION
// ============================================================

/**
 * Get border value for fillMinima algorithm
 */
function borderValue(img) {
  return img.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: roi,
    maxPixels: 10e13
  });
}

/**
 * Improved cloud shadow detection
 * Combines FMask fillMinima with terrain shadow masking
 */
function CloudShadows_S2(img) {
  var nir = img.select('NIR').toInt();
  var swir1 = img.select('SWIR1').toInt();
  
  // FMask fillMinima algorithm
  var finir = ee.Algorithms.FMask.fillMinima(nir, 62480);
  var fiswir1 = ee.Algorithms.FMask.fillMinima(swir1, 40450);
  
  var res1 = finir.subtract(nir);
  var res2 = fiswir1.subtract(swir1);
  
  var result1 = res1.gt(0.02)
    .mask(Water(img).not())
    .mask(PotentialCloud_S2(img).not())
    .unmask(0);
  var result2 = res2.gt(0.02)
    .mask(Water(img).not())
    .mask(PotentialCloud_S2(img).not())
    .unmask(0);
  
  // Potential shadow based on land cover type
  var pshad1 = img.select('bare').eq(60)
    .and(img.select('NIR').lt(0.30))
    .and(img.select('SWIR1').lt(0.28));
  var pshad2 = img.select('bare').neq(60)
    .and(img.select('NIR').lt(0.25))
    .and(img.select('SWIR1').lt(0.18));
  var pshad3 = pshad1.add(pshad2);
  
  // Terrain shadow masking using sun geometry
  var azi = ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE'));
  var zen = ee.Number(img.get('MEAN_SOLAR_ZENITH_ANGLE'));
  var hillshadow = ee.Terrain.hillShadow(img.select('elevation'), azi, zen, 400, true).eq(0);
  
  // Final shadow mask: exclude terrain shadows and cloud pixels
  var shadow1 = pshad3
    .updateMask(hillshadow.not())
    .updateMask(PotentialCloud_S2(img).not())
    .unmask()
    .updateMask(PCPS2(img).not())
    .unmask();
  
  return shadow1.rename('cloudshadow');
}

var shadS2 = S2.map(CloudShadows_S2);
Map.addLayer(shadS2.first(), {palette: ['white', 'blue']}, 'Cloud Shadow');

// ============================================================
// 8. SNOW DETECTION
// ============================================================

/**
 * Multi-threshold snow detection
 * Based on NDSI, NIR reflectance, and Blue/Red ratio
 */
function SnowS2(img) {
  var snow = NDSI(img).gt(0.15)
    .and(img.select('NIR').gt(0.15))
    .and(img.select('Blue').gt(0.28))
    .and(img.select('Blue').divide(img.select('Red')).gt(0.80));
  return snow.rename('snow');
}

var snowS2 = S2.map(SnowS2);
Map.addLayer(snowS2.first(), {palette: ['black', 'white']}, 'Clear-sky Snow');

// ============================================================
// 9. SNOW RECONSTRUCTION UNDER CLOUDS
// ============================================================

/**
 * Elevation-stratified snow reconstruction
 * Estimates snow cover under clouds based on elevation zones
 * 
 * Algorithm:
 * 1. Divide the study area into 300m elevation zones
 * 2. For each zone, calculate mean elevation of clear-sky snow pixels
 * 3. In cloud/shadow covered areas, classify pixels above this elevation as snow
 * 4. Merge with clear-sky snow for complete snow map
 */
function AllsnowALL(img) {
  var snoweleva = SnowS2(img).unmask();
  var a = img.select('elevation').updateMask(snoweleva); // Clear-sky snow elevations
  var allsnow = snoweleva;

  // Combined cloud and shadow mask
  var cloudandcloudshadow = PotentialCloud_S2(img)
    .add(CloudShadows_S2(img))
    .neq(0)
    .unmask();

  // Get elevation range
  var stats = img.select('elevation').reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: roi,
    maxPixels: 10e13
  });

  var max = ee.Number(stats.get('elevation_max'));
  var min = ee.Number(stats.get('elevation_min'));

  // Create elevation zone sequence (300m intervals)
  var steps = ee.List.sequence(min, max.subtract(1), 300);

  // Process each elevation zone
  var snowList = steps.map(function(i) {
    i = ee.Number(i);
    
    // Elevation zone mask
    var aT1 = img.select('elevation').gt(i).and(img.select('elevation').lt(i.add(300)));
    var imgT = img.updateMask(aT1);

    var aT = a.select('elevation').gt(i).and(a.select('elevation').lt(i.add(300)));
    var a1 = img.updateMask(aT);

    // Mean elevation of clear-sky snow in this zone
    var meanheight1 = ee.Number(a1.select('elevation').reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roi,
      maxPixels: 10e13,
      bestEffort: true
    }).get('elevation', 0));

    // Classify cloud-covered pixels above mean snow elevation as snow
    var newsnow = imgT.select('elevation')
      .updateMask(cloudandcloudshadow)
      .gt(meanheight1)
      .unmask()
      .gt(0);
    
    return newsnow;
  });

  // Merge all elevation zones
  var snowImage = ee.ImageCollection.fromImages(snowList).sum().gt(0);
  allsnow = allsnow.add(snowImage).gt(0);

  return allsnow.rename('AllSnow');
}

var allsnowall = S2.map(AllsnowALL);
Map.addLayer(allsnowall.first(), {palette: ['white', 'blue']}, 'Reconstructed Snow (All)');

// ============================================================
// 10. EXPORT RESULTS (Optional)
// ============================================================

// Uncomment to export results to Google Drive
/*
Export.image.toDrive({
  image: TCS2,
  description: 'Cloud_Mask',
  folder: 'SnowRUCCS_Results',
  region: roi,
  scale: 10,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: shadS2.first(),
  description: 'Cloud_Shadow',
  folder: 'SnowRUCCS_Results',
  region: roi,
  scale: 10,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: allsnowall.first(),
  description: 'Snow_Reconstructed',
  folder: 'SnowRUCCS_Results',
  region: roi,
  scale: 10,
  maxPixels: 1e13
});
*/

print('SnowRUCCS processing complete!');
print('Results displayed on map.');
