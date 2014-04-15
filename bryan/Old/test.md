calcLULCByLakeBuffer.R Data Defintions
========================================================
**Last Update=April 14, 2014**
***************
**This data set contains ???????? objects**:

************
FIX:   The following code will open this data definitions file in a browser:
************
  #Data Definitions:
  
    browseURL('https://github.com/jhollist/cyanoLakes/blob/master/bryan/cyanoBioVolData.md')
************
The data can be loaded into R with the following code:
************
  #Get the Data
  
    load(url('https://raw.github.com/jhollist/cyanoLakes/master/bryan/cyanoBioVolData.rda'))
************
**LULC**: raster('L:/Public/Milstead_Lakes/NLCD2006/nlcd_2006_impervious_2011_edition_2014_03_31.img')

raster LULC<-raster('L:/Public/Milstead_Lakes/NLCD2006/nlcd_2006_impervious_2011_edition_2014_03_31.img')

Rows=104,424 Columns=161,190

**Code**  |**Label**  | **Definition**
------------- | ------------- | ------------- 
0|Zero|Not sure what this means yet.  Not in the metadata.  Probably means missing value.
11|Water|Open Water - All areas of open water, generally with less than 25% cover or vegetation or soil
12|IceSnow|Perennial Ice/Snow - All areas characterized by a perennial cover of ice and/or snow, generally greater than 25% of total cover.
21|DevOpen|Developed, Open Space - Includes areas with a mixture of some constructed materials, but mostly vegetation in the form of lawn grasses.  Impervious surfaces account for less than 20 percent of total cover.  These areas most commonly include large-lot single-family housing units, parks, golf courses, and vegetation planted in developed settings for recreation, erosion control, or aesthetic purposes.
22|DevLow|Developed, Low Intensity -Includes areas with a mixture of constructed materials and vegetation.  Impervious surfaces account for 20-49 percent of total cover.  These areas most commonly include single-family housing units.
23|DevMed|Developed, Medium Intensity - Includes areas with a mixture of constructed materials and vegetation. Impervious surfaces account for 50-79 percent of the total cover.  These areas most commonly include single-family housing units.
24|DevHigh|Developed, High Intensity - Includes highly developed areas where people reside or work in high numbers. Examples include apartment complexes, row houses and commercial/industrial.  Impervious surfaces account for 80 to 100 percent of the total cover.
31|Barren|Barren Land (Rock/Sand/Clay) - Barren areas of bedrock, desert pavement, scarps, talus, slides, volcanic material, glacial debris, sand dunes, strip mines, gravel pits and other accumulations of earthen material. Generally, vegetation accounts for less than 15% of total cover.
41|Deciduous|Deciduous Forest - Areas dominated by trees generally greater than 5 meters tall, and greater than 20% of total vegetation cover. More than 75 percent of the tree species shed foliage simultaneously in response to seasonal change.
42|Evergreen|Evergreen Forest - Areas dominated by trees generally greater than 5 meters tall, and greater than 20% of total vegetation cover. More than 75 percent of the tree species maintain their leaves all year. Canopy is never without green foliage.
43|MixedFor|Mixed Forest - Areas dominated by trees generally greater than 5 meters tall, and greater than 20% of total vegetation cover. Neither deciduous nor evergreen species are greater than 75 percent of total tree cover.
52|Shrub|Shrub/Scrub - Areas dominated by shrubs; less than 5 meters tall with shrub canopy typically greater than 20% of total vegetation. This class includes true shrubs, young trees in an early successional stage or trees stunted from environmental conditions.
71|Grass|Grassland/Herbaceous - Areas dominated by grammanoid or herbaceous vegetation, generally greater than 80% of total vegetation.  These areas are not subject to intensive management such as tilling, but can be utilized for grazing.
81|Pasture|Pasture/Hay - Areas of grasses, legumes, or grass-legume mixtures planted for livestock grazing or the production of seed or hay crops, typically on a perennial cycle. Pasture/hay vegetation accounts for greater than 20 percent of total vegetation.
82|Crops|Cultivated Crops - Areas used for the production of annual crops, such as corn, soybeans, vegetables, tobacco, and cotton, and also perennial woody crops such as orchards and vineyards. Crop vegetation accounts for greater than 20 percent of total vegetation. This class also includes all land being actively tilled.
90|WoodyWet|Woody Wetlands - Areas where forest or shrub land vegetation accounts for greater than 20 percent of vegetative cover and the soil or substrate is periodically saturated with or covered with water.
95|HerbWet|Emergent Herbaceous Wetlands - Areas where perennial herbaceous vegetation accounts for greater than 80 percent of vegetative cover and the soil or substrate is periodically saturated with or covered with water.
NA|Ignore|These are pixels from the cropped grid that are outside the mask (buffer) area.  Ignore.

