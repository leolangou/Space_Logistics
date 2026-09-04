# Space Logistics

## Scope & Motivation

This project studies how future **space logistics infrastructure**, particularly refueling assets in geosynchronous orbit (GEO), can be positioned to best support existing satellites.

The goal is to combine orbital data, mathematical structure, and optimization to identify useful logistics locations in GEO.

Project workflow:

https://usarmywestpoint-my.sharepoint.com/:p:/g/personal/leo_langou_westpoint_edu/IQAlXKRWS8_cTKctUxxl5rlsAeXy5EAIVaQiy3-u5RYCvYQ?e=NymldF
## Methodology

```text
Satellite Data
      ↓
Topological Data Analysis
      ↓
Candidate Infrastructure Locations
      ↓
Optimization
      ↓
Testing / Validation
```

The current work uses public GEO satellite position data and TDA methods such as persistent homology and Mapper to characterize the structure of the satellite population. These results will inform optimization models for refueling infrastructure placement and will then be tested against orbital and logistics performance measures.

## File Structure

```text
Space_Logistics/
├── README.md
├── GETTING_STARTED.md
├── requirements.txt
├── us_geo_timeseries.py
├── TDA_test.ipynb
└── data/
```

- `us_geo_timeseries.py` — retrieves and propagates historical U.S. GEO satellite data.
- `TDA_test.ipynb` — analysis and TDA experimentation.
- `requirements.txt` — Python dependencies.
- `GETTING_STARTED.md` — instructions for setting up and running the repository.
- `data/` — generated or locally stored research data.
