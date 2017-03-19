# Measuring transport accessibility


Accessibility is a measure of how easy it is to get places you need to go, including jobs, goods and services, and social and leisure activities. It's also a useful indicator of the effectiveness of a transport system. 

In my research, I compare accessibility by informal transport, such as minibuses, with accessibility by formal public transport, such as bus rapid transit. 


## What this does

`Accessibility.py` is a bit of the code for a fairly simple calculation behind my accessibilty model. Given data on land use and travel time between zones, it computes the accessibility scores for each zone. 

### Input data

Sample of data for small part of Cape Town is in `\sample_data`.

**Land use**

Data on floor area of different land use types for each traffic analysis zone (TAZ) is in `lu_data.csv`. 

**Travel times**

Travel times between zones for different modes: bus rapid transit (called MyCiTi) and minibus taxi in the years 2011 and 2015. Walk times from each zone centroid to the relevant minibus taxi stop are also provided. 

### The accessibility metric

You can choose between two measures of accessibility: 

- **Cumulative opportunities**: the total amount of a given activity type that is reachable within a maximum travel time by a given travel mode. E.g., how many square meters of retail space are reachable by bus within 45 mins.

- **Gravity measure**: the sum of opportunities, positively weighted by attractiveness and negatively weighted by a distance-decaying cost factor. 

(See [Geurs and Van Wee, 2004](http://www.sciencedirect.com/science/article/pii/S0966692303000607) and [El-Geneidy and Levinson, 2006](http://conservancy.umn.edu/handle/11299/638).)

`Accessibility.py` runs the calculation.

## Results

An example of the results mapped for all of Cape Town (not just the sample data) is [here](http://lisarayle.com/ct-accessibility). (UI design is under development - this is just a minimal example.)

