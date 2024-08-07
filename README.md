
The scripts in this folder prepare an exploratory figure for the IPBES NEXUS 
assessment. Specifically they enable the processing and formatting of relevant
scenarios (downloaded in 'raw_data', but not uploaded) and to create (spatial-)temporal visualizations.

# Methodology

The aim of the analysis is to highlight projected trends in several key indicators
that are linked to the various Nexus elements and their "entrypoints" (e.g. Nature conservation)

- Data from ISIMIP, BES-SIM and other sources are compiled and downloaded. 
- A spatial map for the year 2050 was extracted from each source as well as the aggregated
global estimate per time steps. For the aggregation an arithmetic mean was used.
- Indicator values are transformed, if not already done, into relative estimates 
with the baseline year being the year 2015 as starting date.
- Any indicator with inverse scale (e.g. declining intactness) is inverted for the
spatial overlay to emphasize potential linkages.
- All spatial projections for 2050 were normalized (to a scale of 0-1) and then 
averaged per Nexus element (Biodiversity, Food, Water, Health, Climate). 

## Caveats

* This is not a model comparison exercise and it should be noted that there are 
substantial differences among the various projections (as for any model).
* The spatial analysis can only highlight *potential* interlinkages within the 
nexus in the future, not actual ones.
* The data used here are largely from SSP-RCP scenario combinations as used by the IPCC.
Alternative scenario pathways and framing (e.g. Nature Future Framework) were not available
at the time this figure was created.
* Scenarios are not necessarily strictly comparable given that not all indicators
exist per SSP-RCP combination (although see FigureData.xlsx)
* For many nexus elements there might not exist yet relevant projected indicators, 
reflecting a data gap.

## Data sources
See [FigureData](https://github.com/Martin-Jung/IPBES_NexusRelativeScenarioTrends/blob/main/FigureData.xlsx) file.

## License

CC-BY 4.0

## Contact

[Martin Jung](https://iiasa.ac.at/staff/martin-jung)