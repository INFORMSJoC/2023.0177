# Data Descriptions

The instances are taken from the following paper, “Impact of Autonomous Vehicle Assisted Last-Mile Delivery in Urban to Rural Settings, Transportation Science, Vol.56(6), pp.1530-1548 11/2022” by S. Reed, A. M. Campbell, and B. W. Thomas. The instances are posted online and can be downloaded here (https://iro.uiowa.edu/esploro/outputs/dataset/Instances-for-Impact-of-Autonomous-Vehicle/9983903464002771).

For each instance, there are three csv files, which contain the following information.

1. lat: Latitude measurement of customer/depot location
2. lon: Longitude measurement of customer/depot location
3. key: Origin destination pair
4. Distance: Driving (or walking) distance between origin destination pair
6. Time: Time to drive (or walk) between origin destination pair


The following part is taken from the test_instance_readme file from the original repository 

---------------
DATA GENERATION
---------------
These test instances are constructed to reflect changes in customer geography based on their location. The 
United States Department of Agriculture classifies counties on an urban-to-rural continuum (USDA 2013). These test instances 
include one county for each urban-to-rural code in Illinois. In particular, the chosen county has the 
largest population per code. Specifically, the following counties are included:
	1: Cook County (urban)
	2: Winnebago County
	3: Champaign County
	4: La Salle County
	5: Adams County
	6: Fulton County
	7: Jefferson County
	8: Johnson County
	9: Cumberland County (rural).

For each county, we identify areas inside each county that would include service to n = 50 (or 100) customers. 
We take 15% of the county population provided by the USDA as customers (USDA 2013). Assuming the customers are
uniformly distributed over the county, we use the total square miles per county from the National Association
of Counties to determine the population density (NaCo 2010). 

For each county, we generate 10 random square service areas within the county borders defined by
the Illinois Geospatial Clearing House to yield 10 datasets with 50 customers (or 5 datasets with 100 customers 
when applicable) (Illinois 2003). We use the VeRoViz package in Python to generate a uniform distribution of 
locations for the depot and customers within the region (Peng and Murray 2019). These locations are placed on 
the nearest road and defined by their longitude and latitude coordinates. The square service region is divided 
into four quadrants with each customer location being in a single quadrant. The test instances are constructed 
to satisfy the following conditions:
(1) All customer and depot locations are within the given square service area.
(2) There exists at least 10% of the customers in each quadrant of the square service region.

The VeRoViz package uses OpenRouteService as a data provider for generating the walking and driving 
distances/times between locations. OpenRouteService uses data from OpenStreetMap to restrict the driving route 
to the road network and define the pedestrian path for the walking information.


MORE DETAILS
------------
S. Reed, A. M. Campbell, and B. W. Thomas. Impact of Autonomous Vehicle Assisted Last-Mile Delivery in Urban to Rural Settings.

S. Reed, A. M. Campbell, and B. W. Thomas. Does Parking Matter? The Impact of Search Time for Parking on Last-Mile Delivery Optimization


REFERENCES
----------
Illinois 2003. Illinois county boundaries, polygons and lines. Technical report, Illinois Geospatial Data Clearinghouse, 2003.

NACo 2010. County explorer. https://ce.naco.org/.

L. Peng and C. Murray. VeRoViz: A vehicle routing visualization package. https://veroviz.org, 2019.

USDA 2013. Rural-urban continuum codes. Technical report, United States Department of Agriculture (USDA), May 2013.
